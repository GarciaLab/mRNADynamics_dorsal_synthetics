function tfdrivenentryexit_gpu(varargin)

model = "entryexit"; %choices- entryexit, entry, exit, basic
nSims = 1E3; %number of simulations for the kinetic window model
nPlots = 15;
exitOnlyDuringOffStates = true; %determines connectivity of the markov graph
t_cycle = 8;
variableStateNumber = true;
fixKD = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
gpurng(1, "ThreeFry"); %fastest gpu rng

dmax = 4000;
nSilentStates = 1;


[~, dropboxfolder] = getDorsalFolders;
load([dropboxfolder, '\manuscript\window\basic\dataForFitting\DorsalFluoValues.mat'], 'DorsalFluoValues')
% dls = DorsalFluoValues;


dls = logspace(0, log10(dmax), 40);
% kds = logspace(2, 7, nPlots);
% cs = logspace(-5, 2, nPlots);
if fixKD
    kds = 1E3;
else
    kds = logspace(2, 5.5, nPlots);
end
cs = logspace(-1, 3, nPlots);
pi_exits = logspace(-3, 2, nPlots);
pi_entries = logspace(-1, 2, nPlots);
if variableStateNumber
    nentries = [0, 1:2:10];
    moffs = 1:2:10;
else
    nentries = 1;
    moffs = 1;
end
t_cycles = [5, 6, 7, 8, 9, 10];

%%

switch model
    case "entry"
        pi_exits = 0;
    case "basic"
        pi_exits = 0;
        pi_entries = 1E10;
        nentries = 0;
        nSilentStates = 0;
    case "exit"
        cs = logspace(-5, -1, nPlots);
        pi_exits = logspace(-4, 4, nPlots);
        pi_entries = 1E10;
        nentries = 0;
end



params.dls = dls;
params.kds = kds;
params.cs = cs;
params.pi_exits = pi_exits;
params.pi_entries = pi_entries;
params.nentries = nentries;
params.moffs = moffs;
params.t_cycles = t_cycles;
params.exitOnlyDuringOffStates = exitOnlyDuringOffStates;

params.model = model;

if length(params.nentries)==1 && length(params.moffs) == 1
    o = 1;
    p = 1;
    nOffEntryStates = moffs(p) + nentries(o);
    firstoffstate = nentries(o)+1;
    onstate = nentries(o) + moffs(p) +1;
    silentstate = onstate+1;
    nStates = nentries(o) + moffs(p) + 1 + nSilentStates;
    params.nEntryStates = nentries;
    params.nOffStates = moffs;
    params.nStates = nStates;
end


%%
mfpts = gpuArray(nan( length(dls), length(kds), length(pi_exits),...
    length(cs), length(pi_entries), length(nentries),length(moffs),length(t_cycles), 'single'));
factive = gpuArray(zeros( length(dls), length(kds), length(pi_exits),...
    length(cs), length(pi_entries), length(nentries),length(moffs),length(t_cycles), 'single'));
fpts_std = mfpts;

%i'm doing this calculation here to avoid computational overhead in the sim loop
N_cs = length(params.cs);
N_dls = length(params.dls);
N_kds = length(params.kds);
N_pi_exits = length(params.pi_exits);
N_pi_entries = length(params.pi_entries);
N_nentries = length(params.nentries);
N_moffs = length(params.moffs);
N_t_cycles = length(t_cycles);

% progressCounter = 0;
% nParams = numel(factive);
% progress = 0;
% display_interval = round(nParams / 20);

%dls, kds, pi_exits, cs, pi_entries, nentries, moffs
for q = 1:N_t_cycles
    
for o = 1:N_nentries
    
    if nentries(o) ~= 0
        tau_entry = gpuArray(zeros(nentries(o), nSims, numel(pi_entries)));
        
        for k = 1:length(pi_entries)
            tau_entry(:, :, k) = exprnd(pi_entries(k)^-1, [nentries(o), nSims]);
        end
    else
        tau_entry = [];
    end
    
    for p = 1:N_moffs
        
        nOffEntryStates = moffs(p) + nentries(o);
        nStates = nentries(o) + moffs(p) + 1 + nSilentStates;
        
        tau_exit = gpuArray(zeros(nStates-1, nSims, numel(pi_exits), 'single'));
        for k = 1:length(pi_exits)
            tau_exit(:, :, k) = exprnd(pi_exits(k)^-1, [nStates-1, nSims]);
        end
        
        for m = 1:N_cs
            
            
            display("tfdrivenentryexit progress: "+ num2str(( (m-1) / N_cs)*100)+"%" )
            
            %             for i = 1:N_dls
            parfor i = 1:N_dls
                for j = 1:N_kds
                    
                    %pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
                    tau_on = exprnd( ( gpuArray(single(cs(m).* ((dls(i)./kds(j)) ./ (1 + dls(i)./kds(j))))) )^-1,...
                        [moffs(p)+1, nSims]);
                    
                    for n = 1:N_pi_entries
                        
                        if ~isempty(tau_entry)
                            tau_entry_off = [tau_entry(:, :, n); tau_on];
                        else
                            tau_entry_off = tau_on;
                        end
                        
                        for k = 1:N_pi_exits
                            
                            %let's determine if the transition is to the silent
                            %state or other.
                            [~, whichTransition] = min(cat(3,tau_entry_off,tau_exit(:, :, k)), [], 3);
                            
                            if exitOnlyDuringOffStates
                                whichTransition(1:nentries(o), :) = 1;
                            end
                            %the simulations that reached "on" are the ones that
                            %never reached silent.
                            reachedOn = sum(whichTransition(1:nOffEntryStates, :), 1) ==...
                                nOffEntryStates;
                            
                            %let's get the total duration of the successful trajectories up to
                            %the on state
                            onsets_sim = sum(tau_entry_off(1:nOffEntryStates, reachedOn), 1);
                            trunc = onsets_sim < t_cycles(q);
                            
                            factive(i,j,k,m,n, o, p, q) = sum(trunc) / nSims;
                            mfpts(i, j, k, m, n, o, p, q) = mean(onsets_sim(trunc), 'omitnan');
                            %                     fpts_std(i, j, k, m, n, o, p) = std(onsets_sim(trunc));
                            
                        end %for pi_entries
                    end %for pi_exits
                end% for kd
            end %for dl
        end %for c
    end %for moffs
end% for nentries
end % for t_cycles

factive = gather(factive);
mfpts = gather(mfpts);
% 
% dt = mfpts(:, nearestIndex(kds, 1E4), :, :, :, :, :, :) -...
%     mfpts(:, nearestIndex(kds, 400), :, :, :, :, :, :);
% 
% dt = repmat(dt, [1 length(kds) 1 1 1 1 1]);
dt=[];

params.nPoints = 2E4;


saveStr = model;
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
if length(nentries) > 1 || length(moffs) > 1
    saveStr = saveStr + "variableStateNumber";
end

params.saveStr = saveStr;

try
    goodMatrixIndices = plotTFDrivenParams2(factive, dt, mfpts, 'dim', 2, 'params', params);
catch
    goodMatrixIndices = [];
end
%
save(dropboxfolder + "\simulations\" + "tf_paramsearch_"+saveStr+"_.mat",...
    'params', 'mfpts', 'dt', 'factive' , 'goodMatrixIndices')

% %Also save a timestamped copy.
% dttm = strrep(strrep(string(datetime), " ", "_"), ":", "_");
% save(dropboxfolder + "\" + "tf_paramsearch_"+model+"_"+dttm+"_.mat")



% load(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat")

plotTFDrivenParams2(factive, dt, mfpts, 'nPoints', params.nPoints, 'dim', 2, 'params', params)

plotGoodCurves(factive, dt, mfpts, params, goodMatrixIndices)

%%
%
% figure;
% t = tiledlayout('flow');
% for k = 1:2:length(dls)
%     nexttile;
%     try
%         plotTFDrivenParams2(factive(k, :, :, :, :) , dt(k, :, :, :, :),...
%             mfpts(k, :, :, :, :), 'params', params, 'fig', gcf);
%     end
% end
% title(t, 'Effect of [Dorsal] on parameter space');
%
%
% figure;
% t = tiledlayout('flow');
% for k = 1:1:length(cs)
%     nexttile;
%     try
%         plotTFDrivenParams2(factive(:, :, :, k, :) , dt(:, :, :, k, :),...
%             mfpts(:, :, :, k, :),'params', params, 'fig', gcf);
%     end
%     title(t, 'Effect of c on parameter space');
% end
% figure;
% t = tiledlayout('flow');
% for k = 1:1:length(kds)
%     nexttile;
%     try
%         plotTFDrivenParams2(factive(:, k, :, :, :) ,...
%             dt(:, k, :, :, :), mfpts(:, k, :, :, :),'params', params,'nPoints', 1E3, 'fig', gcf);
%     end
% end
% title(t, 'Effect of KD on parameter space');
%
% figure;
% t = tiledlayout('flow');
% for k = 1:1:length(pi_exits)
%     nexttile;
%     try
%         plotTFDrivenParams2(factive(:, :, k, :, :) ,...
%             dt(:, :, k, :, :), mfpts(:, :, k, :, :),'params', params, 'nPoints', 1E3, 'fig', gcf);
%         %          title(['\pi_{exit} = ', num2str(round2(pi_exits(k))), ' min^{-1}'])
%         title(num2str(round2(pi_exits(k))))
%     end
% end
% title(t, 'Effect of pi_exit on parameter space');
%
%
%
% figure;
% t = tiledlayout('flow');
% for k = 1:1:length(pi_entries)
%     nexttile;
%     try
%         plotTFDrivenParams2(factive(:, :, :, :, k),...
%             dt(:, :, :, :, k), mfpts(:, :, :, :, k), 'params', params, 'fig', gcf,  'nPoints', 1E3);
%     end
%     %         title(['\pi_{entry} = ', num2str(round2(pi_entries(k))), ' min^{-1}'])
%     title(num2str(round2(pi_entries(k))))
% end
% title(t, 'Effect of pi_entry on parameter space');
%
%


figure;
t = tiledlayout('flow');
for k = 1:3:length(params.pi_exits)
    nexttile;
    try
        plotTFDrivenParams2(factive(:, :, :, :, k, :, :),...
            dt(:, :, :, :, k, :, :), mfpts(:, :, :, :, k, :, :), 'params', params, 'fig', gcf,  'nPoints', 1E3);
    end
    %         title(['\pi_{entry} = ', num2str(round2(pi_entries(k))), ' min^{-1}'])
    title(num2str(round2(params.pi_exits(k))))
end
title(t, 'Effect of pi_entry on parameter space');


disp('debug stop')
end

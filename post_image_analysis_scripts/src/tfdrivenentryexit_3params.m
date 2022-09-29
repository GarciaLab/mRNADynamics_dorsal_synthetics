%leaving this here for to remember for testing
%setpref('profiler','showJitLines',1);
%profile -memory on;

tic

% model = "basic";
% model = "entry";
% model = "exit";
model = "entryexit";

rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
dmax = 4000;


t_cycle = 8; %min

nOffStates = 5;
nEntryStates = 5;
nSilentStates = 1;

occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );

% %% SLOW
% exitOnlyDuringOffStates = true;
% nSims = 1E5;
% nPlots = 20;
% 
% dls = logspace(0, log10(dmax), 80);
% kds = logspace(2, 7, nPlots);
% cs = logspace(-5, 2, nPlots);
% pi1s = logspace(-3, 2, nPlots);
% pi2s = logspace(-1, 2, nPlots);

% %% MEDIUM
% exitOnlyDuringOffStates = true;
% nSims = 1E4;
% nPlots = 20;
% 
% dls = logspace(0, log10(dmax), 20);
% kds = logspace(2, 7, nPlots);
% cs = logspace(-5, 2, nPlots);
% pi1s = logspace(-3, 2, nPlots);
% pi2s = logspace(-1, 2, nPlots);
% 
% %%

%% FAST
exitOnlyDuringOffStates = true;
nSims = 2
nPlots = 3;

pi_basics = logspace(-3, 2, nPlots);
pi_exits = logspace(-3, 2, nPlots);
pi_entries = logspace(-1, 2, nPlots);
%%

switch model
    case "entry"
        pi_exits = 0;
    case "basic"
        pi_exits = 0;
        pi_entries = 1E10;
    case "exit"
        pi_exits = logspace(-4, 4, nPlots);
        pi_entries = 1E10;
end

nOffEntryStates = nOffStates + nEntryStates;
firstoffstate = nEntryStates+1;
onstate = nEntryStates + nOffStates+1;
silentstate = onstate+1;
nStates = nEntryStates + nOffStates + 1 + nSilentStates;

% nParams = numel(dls)*numel(kds)*numel(pi1s)*numel(pi2s)*numel(cs);

params.model = model;
params.nEntryStates = nEntryStates;
params.nOffStates = nOffStates;
params.nStates = nStates;
params.exitOnlyDuringOffStates = exitOnlyDuringOffStates;
params.pi1s = pi_exits;
params.pi2s = pi_entries;
params.pibasics = pi_basics;
params.model = model;
params.nEntryStates = nEntryStates;
params.nOffStates = nOffStates;
params.nStates = nStates;
params.exitOnlyDuringOffStates = exitOnlyDuringOffStates;


%%
mfpts = nan(length(pi_basics), length(pi_exits), length(pi_entries));
factive =  nan(length(pi_basics), length(pi_exits), length(pi_entries));
fpts_std = nan(length(pi_basics), length(pi_exits), length(pi_entries));

N_pi1s = length(params.pi1s);
N_pi2s = length(params.pi2s);
N_pibasics = length(params.pibasics); 

tau_on = nan(nEntryStates+1, nSims, numel(pi_basics), 'double');
for k = 1:length(pi_basics)
    tau_on(:, :, k) = exprnd(pi_basics(k)^-1, [nOffStates+1, nSims]);
end

tau_exit = nan(nStates-1, nSims, numel(pi_exits), 'double');
for k = 1:length(pi_exits)
    tau_exit(:, :, k) = exprnd(pi_exits(k)^-1, [nStates-1, nSims]);
end

tau_entry = nan(nEntryStates, nSims, numel(pi_entries), 'double');
for k = 1:length(pi_entries)
    tau_entry(:, :, k) = exprnd(pi_entries(k)^-1, [nEntryStates, nSims]);
end

tau_entry_off = cat(1, tau_entry, tau_on); 

 [~, whichTransition] = min(cat(4,tau_entry_off, tau_exit), [], 4);
 
if exitOnlyDuringOffStates
    whichTransition(1:nEntryStates, :, :) = 1;
end

onsets_sim = squeeze(sum(tau_entry_off, 1));

reachedOn = squeeze(sum(whichTransition(1:nOffEntryStates, :, :), 1) ==...
                        nOffEntryStates);     
 trunc = onsets_sim < t_cycle;
                   
                    
 a = tau_entry_off;
 b = find(~reachedOn)
 a(:, b) =[];
tau_entry_off_reachedOn = tau_entry_off(1:nOffEntryStates, reachedOn);
onsets_sim = sum(tau_entry_off_reachedOn, 1);

onsets_sim = sum(tau_entry_off(1:nOffEntryStates, reachedOn), 1);
 
                                      
factive = sum(reachedOn(trunc), 1)/nSims;

mfpts(i, k, n) = mean(onsets_sim(trunc), 1);
% display("tfdrivenentryexit progress: "+ num2str(( (m-1) / N_cs)*100)+"%" )

%pi_basics, pi_exits,pi_entries
for i = 1:N_pibasics

            for n = 1:N_pi2s
                                
                for k = 1:N_pi1s
                    
                    %let's determine if the transition is to the silent
                    %state or other.
                    [~, whichTransition] = min(cat(3,tau_entry_off,tau_exit(:, :, k)), [], 3);
                    
                    if exitOnlyDuringOffStates
                        whichTransition(1:nEntryStates, :) = 1;
                    end
                    %the simulations that reached "on" are the ones that
                    %never reached silent.
                    reachedOn = sum(whichTransition(1:nOffEntryStates, :), 1) ==...
                        nOffEntryStates;
                    
                    %let's get the total duration of the successful trajectories up to
                    %the on state
                    onsets_sim = sum(tau_entry_off(1:nOffEntryStates, reachedOn), 1);
                    trunc = onsets_sim < t_cycle;
                                      
                    factive(i,k, n) = sum(reachedOn(trunc))/nSims;
                    mfpts(i, k, n) = mean(onsets_sim(trunc));
                    %                     fpts_std(i, j, k, m, n) = std(onsets_sim(trunc));
                    
                end %for pi2
            end %for pi1
        end% for pibasic

dt = mfpts(:, nearestIndex(kds, 1E4), :, :, :) -...
    mfpts(:, nearestIndex(kds, 400), :, :, :);

dt = repmat(dt, [1 length(kds) 1 1 1]);


params.nPoints = 2E4; 

[~, dropboxfolder] = getDorsalFolders;

saveStr = model + "_3param_";
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
params.saveStr = saveStr;

goodMatrixIndices = plotTFDrivenParams2(factive, dt, mfpts, 'dim', 2, 'params', params);

save(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat",...
    'params', 'mfpts', 'dt', 'factive' , 'goodMatrixIndices')

% %Also save a timestamped copy.
% dttm = strrep(strrep(string(datetime), " ", "_"), ":", "_");
% save(dropboxfolder + "\" + "tf_paramsearch_"+model+"_"+dttm+"_.mat")


toc
% load(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat")

plotTFDrivenParams2(factive, dt, mfpts, 'nPoints', nPoints, 'dim', 2, 'params', params)

plotGoodCurves(factive, dt, mfpts, params, goodMatrixIndices)

%%

    figure;
    t = tiledlayout('flow');
    for k = 1:2:length(dls)
        nexttile;
        try
        plotTFDrivenParams2(factive(k, :, :, :, :) , dt(k, :, :, :, :),...
            mfpts(k, :, :, :, :), 'params', params, 'fig', gcf);
        end
    end
    title(t, 'Effect of [Dorsal] on parameter space');


    figure;
    t = tiledlayout('flow');
    for k = 1:1:length(cs)
        nexttile;
        try
        plotTFDrivenParams2(factive(:, :, :, k, :) , dt(:, :, :, k, :),...
            mfpts(:, :, :, k, :),'params', params, 'fig', gcf);
        end
    title(t, 'Effect of c on parameter space');
    end
    figure;
    t = tiledlayout('flow');
    for k = 1:1:length(kds)
        nexttile;
        try
            plotTFDrivenParams2(factive(:, k, :, :, :) ,...
                dt(:, k, :, :, :), mfpts(:, k, :, :, :),'params', params,'nPoints', 1E3, 'fig', gcf);
        end
    end
    title(t, 'Effect of KD on parameter space');

    figure;
    t = tiledlayout('flow');
    for k = 1:1:length(pi1s)
        nexttile;
        try
            plotTFDrivenParams2(factive(:, :, k, :, :) ,...
                dt(:, :, k, :, :), mfpts(:, :, k, :, :),'params', params, 'nPoints', 1E3, 'fig', gcf);
            %          title(['\pi_{exit} = ', num2str(round2(pi1s(k))), ' min^{-1}'])
            title(num2str(round2(pi1s(k))))
        end
    end
    title(t, 'Effect of pi_exit on parameter space');



    figure;
    t = tiledlayout('flow');
    for k = 1:1:length(pi2s)
        nexttile;
            try
            plotTFDrivenParams2(factive(:, :, :, :, k),...
                dt(:, :, :, :, k), mfpts(:, :, :, :, k), 'params', params, 'fig', gcf,  'nPoints', 1E3);
            end
        %         title(['\pi_{entry} = ', num2str(round2(pi2s(k))), ' min^{-1}'])
        title(num2str(round2(pi2s(k))))
    end
    title(t, 'Effect of pi_entry on parameter space');

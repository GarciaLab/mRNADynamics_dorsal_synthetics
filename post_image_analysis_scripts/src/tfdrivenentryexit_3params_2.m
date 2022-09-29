function tfdrivenentryexit_3params_2(varargin)

% model = "basic";
% model = "entry";
% model = "exit";
model = "entryexit";

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

%leaving this here for to remember for testing
%setpref('profiler','showJitLines',1);
%profile -memory on;

tic


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
nSims = 1E3;
nPlots = 200;

% pi_basics = logspace(-1, 1, nPlots);
pi_basics = logspace(-3, 15, nPlots);
pi_exits = logspace(-3, 2, nPlots);
% pi_exits = 0;
% pi_entries = logspace(-1, 1, nPlots);
pi_entries = logspace(-3, 15, nPlots);

% pi_entries = 1E10;
%%

switch model
    case "entry"
        pi_exits = 0;
    case "basic"
        pi_exits = 0;
        pi_entries = 1E10;
    case "exit"
%         pi_exits = logspace(-4, 4, nPlots);
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
params.pi_exits = pi_exits;
params.pi_entries = pi_entries;
params.pi_basics = pi_basics;
params.model = model;
params.nEntryStates = nEntryStates;
params.nOffStates = nOffStates;
params.nStates = nStates;
params.exitOnlyDuringOffStates = exitOnlyDuringOffStates;


%%
mfpts = nan(length(pi_basics), length(pi_entries), length(pi_exits));
factive =  nan(length(pi_basics), length(pi_entries), length(pi_exits));
fpts_std = nan(length(pi_basics), length(pi_entries), length(pi_exits));

N_pi_exits = length(params.pi_exits);
N_pi_entries = length(params.pi_entries);
N_pi_basics = length(params.pi_basics); 

tau_on = nan(nOffStates+1, nSims, numel(params.pi_basics), 'double');
for k = 1:length(pi_basics)
    tau_on(:, :, k) = exprnd(pi_basics(k)^-1, [nOffStates+1, nSims]);
end

tau_exit = nan(nStates-1, nSims, numel(params.pi_exits), 'double');
for k = 1:length(pi_exits)
    tau_exit(:, :, k) = exprnd(pi_exits(k)^-1, [nStates-1, nSims]);
end

tau_entry = nan(nEntryStates, nSims, numel(params.pi_entries), 'double');
for k = 1:length(pi_entries)
    tau_entry(:, :, k) = exprnd(pi_entries(k)^-1, [nEntryStates, nSims]);
end

% if contains(model, "entry")
%     tau_entry_off = cat(1, tau_entry, tau_on); 
% else
%     tau_entry_off = tau_on;
% end
% 


%pi_basics, pi_entries, pi_exits
for i = 1:N_pi_basics
    
display("tfdrivenentryexit progress: "+ num2str(( (i-1) / N_pi_basics)*100)+"%" )

            for j = 1:N_pi_entries
                        
                tau_entry_off = [tau_entry(:, :, j); tau_on(:, :, i)];
                
                for k = 1:N_pi_exits
                    
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
                                      
                    factive(i, j, k) = sum(reachedOn(trunc))/nSims;
                    mfpts(i, j, k) = mean(onsets_sim(trunc));
                    %                     fpts_std(i, j, k, m, n) = std(onsets_sim(trunc));
                    
                end %for pi_exit
            end %for pi_entry
        end% for pibasic

params.nPoints = 2E4; 

[~, dropboxfolder] = getDorsalFolders;

saveStr = model + "_3param_";
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
params.saveStr = saveStr;

goodMatrixIndices = plotTFDrivenParams_3params(factive, mfpts, 'params', params);

save(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat",...
    'params', 'mfpts',  'factive' , 'goodMatrixIndices')

% %Also save a timestamped copy.
% dttm = strrep(strrep(string(datetime), " ", "_"), ":", "_");
% save(dropboxfolder + "\" + "tf_paramsearch_"+model+"_"+dttm+"_.mat")


toc
% load(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat")

plotTFDrivenParams_3params(factive, mfpts, 'nPoints', params.nPoints,  'params', params)

% plotGoodCurves(factive,mfpts, params, goodMatrixIndices)

% %%
% try
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:5:length(params.pibasics)
%         nexttile;
% %         try
%             plotTFDrivenParams_3params(factive(k, :, :) ,...
%                 mfpts(k, :, :),'params', params, 'nPoints', 1E3, 'fig', gcf);
%             %          title(['\pi_{exit} = ', num2str(round2(pi1s(k))), ' min^{-1}'])
%             title(num2str(round2(pi_basics(k))))
% %         end
%     end
%     title(t, 'Effect of pi_basic on parameter space');
%     
%     
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:5:length(pi_exits)
%         nexttile;
%         try
%             plotTFDrivenParams_3params(factive(:, :, k) ,...
%                 mfpts(:, :, k),'params', params, 'nPoints', 1E3, 'fig', gcf);
%             %          title(['\pi_{exit} = ', num2str(round2(pi1s(k))), ' min^{-1}'])
%             title(num2str(round2(pi_exits(k))))
%         end
%     end
%     title(t, 'Effect of pi_exit on parameter space');
% 
% 
% 
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:5:length(pi_entries)
%         nexttile;
%             plotTFDrivenParams_3params(factive(:,k,:),...
%                 mfpts(:,k,:), 'params', params, 'fig', gcf,  'nPoints', 1E3);
%         %         title(['\pi_{entry} = ', num2str(round2(pi2s(k))), ' min^{-1}'])
%         title(num2str(round2(pi_entries(k))))
%     end
%     title(t, 'Effect of pi_entry on parameter space');
    
%     
%       figure;
%     for k = 1:2:length(params.pi2s)
%             plotTFDrivenParams_3params(factive(:,k,nearestIndex(params.pi1s, 1)),...
%                 mfpts(:,k,nearestIndex(params.pi1s, 1)), 'params', params, 'fig', gcf,  'nPoints', 1E3);
%         %         title(['\pi_{entry} = ', num2str(round2(pi2s(k))), ' min^{-1}'])
% %         title(num2str(round2(pi_entries(k))))
%         hold on
%     end
% %     title(t, 'Effect of pi_entry on parameter space');
% % end


figure;
tiledlayout('flow')
for j = 1:15:length(params.pi_exits)
    nexttile;
    for k = 1:15:length(params.pi_entries)
            plotTFDrivenParams_3params(factive(:,k,j),...
                mfpts(:,k,j), 'params', params, 'fig', gcf,  'nPoints', 1E3);
        %         title(['\pi_{entry} = ', num2str(round2(pi2s(k))), ' min^{-1}'])
%         title(num2str(round2(params.pi_entries(k))))
        hold on
    end
    title(num2str(round2(params.pi_exits(j))))
end
% %     
%     figure;
%     for k = 1:2:length(pi_entries)
%             plotTFDrivenParams_3params(factive(:,k,1),...
%                 mfpts(:,k,1), 'params', params, 'fig', gcf,  'nPoints', 1E3);
%         %         title(['\pi_{entry} = ', num2str(round2(pi2s(k))), ' min^{-1}'])
%         title(num2str(round2(pi_entries(k))))
%         hold on
%     end
%     title(t, 'Effect of pi_entry on parameter space');

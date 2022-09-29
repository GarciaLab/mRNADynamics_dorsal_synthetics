function generateMarkovTrajectoriesWrapper(model)
% 
% model = "entry";
% model = "entryexit"; 
%model = "basic"

rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
dmax = 4000;


t_cycle = 8; %min

nOffStates = 5;
nEntryStates = 5;
nSilentStates = 1;

occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );


exitOnlyDuringOffStates = true;
nSims = 30;
nPlots = 1;



dls = 1000;
kds = 1E3;
cs = 10;
pi1s = .5;
pi2s = 1;

switch model
    case "entry"
        pi1s = 0;
    case "basic"
        cs = 1;
        pi1s = 0;
        pi2s = 1E10;
        nEntryStates = 0;
    case "exit"
        pi2s = 1E10;
end

nOffEntryStates = nOffStates + nEntryStates;
firstoffstate = nEntryStates+1;
onstate = nEntryStates + nOffStates+1;
silentstate = onstate+1;
nStates = nEntryStates + nOffStates + 1 + nSilentStates;

tau_exit = nan(nStates-1, nSims, numel(pi1s), 'double');
for k = 1:length(pi1s)
    tau_exit(:, :, k) = exprnd(pi1s(k)^-1, [nStates-1, nSims]);
end

tau_entry = nan(nEntryStates, nSims, numel(pi2s), 'double');
for k = 1:length(pi2s)
    tau_entry(:, :, k) = exprnd(pi2s(k)^-1, [nEntryStates, nSims]);
end

clear params;
params.dls = dls;
params.kds = kds;
params.cs = cs;
params.pi1s = pi1s;
params.pi2s = pi2s;
params.model = model;
params.nEntryStates = nEntryStates;
params.nOffStates = nOffStates;
params.nStates = nStates;
params.exitOnlyDuringOffStates = exitOnlyDuringOffStates;


%%
mfpts = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));
factive = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));
fpts_std = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));

N_cs = length(params.cs);
N_dls = length(params.dls);
N_kds = length(params.kds);
N_pi1s = length(params.pi1s);
N_pi2s = length(params.pi2s);

%dls, kds, pi1s, cs, pi2s
for m = 1:N_cs
    
    
    for i = 1:N_dls
        for j = 1:N_kds
            
            pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
            tau_on = exprnd(pi0^-1, [nOffStates+1, nSims]);
            
            for n = 1:N_pi2s
                
                tau_entry_off = [squeeze(tau_entry(:, :, n)); tau_on];
                
                for k = 1:N_pi1s
                    
                    %let's determine if the transition is to the silent
                    %state or other.
                    [~, whichTransition] = min(cat(3,tau_entry_off,squeeze(tau_exit(:, :, k))), [], 3);
                    
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
                    onsets_sim_truncated = onsets_sim(onsets_sim < t_cycle);
                                     
                    
                end %for pi2
            end %for pi1
        end% for kd
    end %for dl
end %for c

generateMarkovTrajectories(whichTransition, tau_entry_off, tau_exit, model);
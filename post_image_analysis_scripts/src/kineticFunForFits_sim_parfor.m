function y = kineticFunForFits_sim_parfor(x, theta, modelOpts)

%sims- dls, kds, pi1s (piexit), cs, pi2s(pientry), nentries, moffs
%theta- c, kd, n, m, pientry, piexit
%example input:
%y = kineticFunForFits_sim(10:250:4000,  [2.5E5, 3300, 5, 5, 13.5, .87], []);

rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period

if isempty(modelOpts)
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = 1E3;
end


nSims = modelOpts.nSims;
exitOnlyDuringOffStates = modelOpts.exitOnlyDuringOffStates;
t_cycle = 8; %min

%mcmcpred feeds in x as a struct. otherwise it's an array
if isstruct(x)
    x = x.ydata(:, 1);
end

dls = x;
kds = theta(2);
cs = theta(1);
nentries = round(theta(3));
moffs = round(theta(4));
pi_entries = theta(5);
pi_exits = theta(6);

nSilentStates = 1;
nOffEntryStates = moffs + nentries;
nStates = nentries + moffs+ 1 + nSilentStates;

if nentries > 0 
    tau_entry = exprnd(gpuArray(pi_entries^-1), [nentries, nSims]);
else
    tau_entry = [];
end

occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );
% occupancy2 = @(d,kd) d.*kd

tau_on = zeros(moffs + 1, nSims, length(dls), 'double');
parfor k = 1:length(dls)
    tau_on(:, :, k) = exprnd(  ((cs.*occupancy(dls(k), kds)).^-1), [moffs+1, nSims]);
end

tau_exit = exprnd(gpuArray(pi_exits^-1), [nStates-1, nSims]);

factive = nan(length(dls), 1);
onset = factive;
for k = 1:length(dls)
    
    tau_entry_off = [tau_entry; tau_on(:, :, k)];
    
    %let's determine if the transition is to the silent
    %state or other.
    if nSilentStates == 1
        [~, whichTransition] = min(cat(3,tau_entry_off,tau_exit), [], 3);
        if exitOnlyDuringOffStates
            whichTransition(1:nentries, :) = 1;
        end
    else
        whichTransition = ones(size(tau_entry_off));
    end
    
    
    %the simulations that reached "on" are the ones that
    %never reached silent.
    reachedOn = sum(whichTransition(1:nOffEntryStates, :), 1) ==...
        nOffEntryStates;
    
    %let's get the total duration of the successful trajectories up to
    %the on state
    onsets_sim = sum(tau_entry_off(1:nOffEntryStates, reachedOn), 1);
    trunc = onsets_sim < t_cycle;
    
    factive(k) = sum(reachedOn(trunc))/nSims;
    onset(k) =  mean(onsets_sim(trunc));
end

y = [factive, onset];

end
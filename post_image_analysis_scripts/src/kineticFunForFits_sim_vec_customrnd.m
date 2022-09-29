function y = kineticFunForFits_sim_vec_customrnd(x, theta, modelOpts)

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
n_dls = length(dls);
kds = theta(2);
cs = theta(1);
nentries = theta(3);
moffs = theta(4);
pi_entries = theta(5);
pi_exits = theta(6);

nSilentStates = 1;
nOffEntryStates = moffs + nentries;
nStates = nentries + moffs+ 1 + nSilentStates;

%generate enough random numbers up front for the whole sim
r_vec = rand(1, nentries*nSims + ( (moffs+1) * nSims * n_dls ) + ( (nStates-1) * nSims));

mu = pi_entries^-1;
sizeOut = [nentries, nSims];
n1 = sizeOut(1) * sizeOut(2); 
r_mat = reshape(r_vec( 1: sizeOut(1)*sizeOut(2) ), sizeOut);
tau_entry = repmat(-mu .* log(r_mat), 1, 1, n_dls);

% tau_entry = exprnd( (pi_entries^-1), [nentries, nSims]);

mu = pi_exits^-1;
sizeOut = [nStates-1, nSims];
n2 = sizeOut(1) * sizeOut(2); 
r_mat = reshape(r_vec( n1 + 1: n1+ n2 ), sizeOut);
tau_exit = repmat(-mu .* log(r_mat), 1, 1, n_dls);


pi_on_inv = ( (cs .* ( (dls./kds) ./ (1 + dls./kds) ) ).^-1)';
tau_on = zeros(moffs + 1, nSims, n_dls, 'double');
for k = 1:n_dls   
    mu =  pi_on_inv(k);
    sizeOut = [moffs+1, nSims];
    n3 = sizeOut(1) * sizeOut(2); 
    r_mat = reshape(r_vec( n1 + n2 + 1: n1 + n2+ n3 ), sizeOut);
    tau_on(:, :, k) = -mu .* log(r_mat);
end


tau_entry_off = cat(1, tau_entry, tau_on);

%let's determine if the transition is to the silent
%state or other.
[~, whichTransition] = min(cat(4,tau_entry_off,tau_exit), [], 4);
if exitOnlyDuringOffStates
    whichTransition(1:nentries, :, :) = 1;
end

%the simulations that reached "on" are the ones that
%never reached silent.
reachedOn = sum(whichTransition(1:nOffEntryStates, :, :), 1) ==...
        nOffEntryStates;

%let's get the total duration of the successful trajectories up to
%the on state
onsets_sim= squeeze(sum(tau_entry_off(1:nOffEntryStates, :, :) .*  repmat(reachedOn, nOffEntryStates, 1, 1), 1));
onsets_sim(onsets_sim == 0) = nan; 
trunc = onsets_sim < t_cycle;
 
%finally, we'll compute the fraction active and mean first passage time to
%the on state. 
factive = sum(trunc, 1) /nSims;

temp = (onsets_sim .* trunc);
temp(temp==0) = nan;
onset =  mean(temp, 1, 'omitnan');
    
y = [factive', onset'];

end
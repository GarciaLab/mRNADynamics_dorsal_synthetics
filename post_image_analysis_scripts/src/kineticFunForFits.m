function y = kineticFunForFits_table(x, theta, modelOpts)

sims = modelOpts.sims; 

%sims- dls, kds, pi1s (piexit), cs, pi2s(pientry)
%theta- c, kd, n, m, pientry, piexit

%mcmcpred feeds in x as a struct. otherwise it's an array
if isstruct(x)
    x = x.ydata(:, 1);
end

for k = 1:length(x)
    ind_dl(k) = nearestIndex(sims.params.dls, x(k));
end
ind_kd = nearestIndex(sims.params.kds ,theta(2) );
ind_c = nearestIndex(sims.params.cs, theta(1));
ind_pientry = nearestIndex(sims.params.pi1s, theta(5));
ind_piexit = nearestIndex(sims.params.pi2s, theta(6));

factive = sims.factive(ind_dl, ind_kd, ind_pientry, ind_c, ind_piexit); 
onset = sims.mfpts(ind_dl, ind_kd, ind_pientry, ind_c, ind_piexit); 
% 
% if any(factive > .5)
%     'stop'
% end
y = [factive, onset];

end
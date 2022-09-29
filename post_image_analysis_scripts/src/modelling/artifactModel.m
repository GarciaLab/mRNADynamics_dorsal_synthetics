close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])
a = combinedCompiledProjects_allEnhancers(...
    [combinedCompiledProjects_allEnhancers.cycle]==12);


figure;
b = a;
fluos = [];
for k = 1:length(b)
    fluos = [fluos, max(b(k).particleFluo3Slice)];
end  
% fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
offset = min(fluos) - 1;
fluos = fluos - offset;
R_max = prctile(fluos, 95);
fluos = fluos./R_max;
fluos = log(fluos);
 if ~isempty(fluos)
    histogram(fluos, 'Normalization', 'pdf', 'facealpha', .5)
    hold on;
%     pd = fitdist(fluos','Lognormal');
    pd = fitdist(fluos','Normal');
%     x_values = .01:.01:3000;
    x_values = min(fluos)*.9:.01:max(fluos);
    y = pdf(pd,x_values);
    plot(x_values,y,'-','LineWidth',1.5)
 end
 coeffText = getDistributionText(pd);
 
% title(['max spot fluorescence. Lognormal fit';coeffText'])
title(['normalized max spot fluorescence. Normal fit';coeffText'])
ylabel('pdf')
xlabel(' log normalized max fluo (au)')

%%

% v = .67;
% v = .41;
% R = prctile(fluos, 95);
R = 1;
% v = pd.sigma;
% v = .41
Lo = (54 - offset) ./R_max;
Ls  = R.*([.05, .1, .25, .9]);
% L = R.*.25;
L = Lo;
kd = 810;
kds = [250, 500, 1000, 5000];
d = [0.01:1:4000];
n = 1;

dmRNAdt = @(d, kd, n, R) (1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R;

dmRNAdtDist = @(d, kd, n, v, R, L) (0.1E1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R.*( ...
  0.1E1+erf(0.353553E0.*v.^(-1).*(v.^2+(-0.2E1).*log(L)+0.2E1.*log( ...
  R+(-0.1E1).*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R)))).*erfc( ...
  0.353553E0.*v.^(-1).*(v.^2+0.2E1.*log(L)+(-0.2E1).*log(R+(-0.1E1) ...
  .*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R))).^(-1);

dmRNAdtObs = @(d, kd, n, v, R, L) (0.1E1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R.*( ...
  0.1E1+erf(0.353553E0.*v.^(-1).*(v.^2+(-0.2E1).*log(L)+0.2E1.*log( ...
  R+(-0.1E1).*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R)))).*erfc( ...
  0.353553E0.*v.^(-1).*(v.^2+0.2E1.*log(L)+(-0.2E1).*log(R+(-0.1E1) ...
 .*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R))).^(-1);

fractionActive = @(d, kd, n, v, R, L)  (1/2).*erfc((-1).*2.^(-1/2).*v.^(-1).*((-0.5E0).*v.^2+( ...
  -1).*log(L)+log((1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R)));


figure; tiledlayout('flow')
for k = 1:length(Ls)
    nexttile;
    yyaxis left
    plot(d, dmRNAdt(d, kd, n, R),'LineWidth', 2 );
    hold on
    plot(d, dmRNAdtObs(d, kd, n, v, R, Ls(k)), '-g', 'LineWidth', 2)
    yyaxis right
    plot(d, fractionActive(d, kd, n, v, R, Ls(k)), 'LineWidth', 2);
    xlabel('[Dorsal] (au)')
    legend({'dmRNA/dt true', 'dmRNA/dt observed', 'fraction active'});
    title("L = " + Ls(k) );
end


figure; tiledlayout('flow')
for k = 1:length(kds)
    nexttile;
    yyaxis left
    plot(d, dmRNAdt(d, kds(k), n, R),'LineWidth', 2 );
    hold on
    plot(d, dmRNAdtObs(d, kds(k), n, v, R, L),'-g','LineWidth', 2 );
    yyaxis right
    plot(d, fractionActive(d, kds(k), n, v, R, L), 'LineWidth', 2);
    xlabel('[Dorsal] (au)')
    legend({'dmRNA/dt true', 'dmRNA/dt observed', 'fraction active'});
    title("KD = " + kds(k) );
end


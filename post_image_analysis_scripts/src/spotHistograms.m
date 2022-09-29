close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])
% % 
% a = combinedCompiledProjects_allEnhancers(strcmpi({combinedCompiledProjects_allEnhancers.dataSet}, '1Dg11_2xDl' )  &...
%     [combinedCompiledProjects_allEnhancers.cycle]==12);

% a = combinedCompiledProjects_allEnhancers(strcmpi({combinedCompiledProjects_allEnhancers.dataSet}, '1Dg11_2xDl' )  &...
%     [combinedCompiledProjects_allEnhancers.cycle]==12 &  [combinedCompiledProjects_allEnhancers.dorsalFluoBin]==11);
% 
a = combinedCompiledProjects_allEnhancers(...
    [combinedCompiledProjects_allEnhancers.cycle]==12);

% b = a(cellfun(@any, {a.particleFrames}))

%%
try
fig = figure;
ax = axes(fig);
maxBin = max([a.dorsalFluoBin]);
for bin = 1:maxBin
    b = a([a.dorsalFluoBin]==bin);
    fluos = [b.particleFluo3Slice];
    fluos(fluos <= 0) = [];
    fluos(isnan(fluos)) = [];
 if ~isempty(fluos)
    histogram(log(fluos), 'Normalization', 'pdf', 'facealpha', .5)
    pd = fitdist(log(fluos)','Normal');
    x_values = 0:.1:10;
    y = pdf(pd,x_values);
    plot(ax, x_values,y,'-','LineWidth',.5)
    hold on
 end
end
title('log spot fluorescence divided by Dl bin. Normal fit')
ylabel('pdf')
end


%%
RNAPFluo = 8.837; %fluorescence AUs of a single RNAP according to our P2P calibration

fig2 = figure;
ax2 = axes(fig2);
b = a;
fluos = [];

for k = 1:length(b)
    fluos = [fluos, min(b(k).particleFluo3Slice)];
end  

fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
if ~isempty(fluos)
    histogram(ax2, fluos, 'Normalization', 'pdf', 'facealpha', .5)
    hold on;
    %   pd = fitdist(fluos','Lognormal');
      pd = fitdist(fluos','Gamma');
    x_values = .1:.1:400;
    y = pdf(pd,x_values);
    plot(ax2, x_values,y,'-','LineWidth',.5)
    plot([mean(fluos) mean(fluos)],[0 0.015],'k-')
end
coeffText = getDistributionText(pd);
title(['min spot fluorescence. Lognormal fit. Mean = ' num2str(mean(fluos));coeffText'])
ylabel('pdf')
xlim([-10 350])
hold off

figure
hold on
AbsFluos = fluos./RNAPFluo;
histogram(AbsFluos, 'Normalization', 'pdf', 'facealpha', .5)
plot([mean(AbsFluos) mean(AbsFluos)],[0 0.15],'k-')
xlim([-10 350]./RNAPFluo)
hold off
title(['Mean = ' num2str(mean(AbsFluos))])



%%

fig3 = figure;
ax3 = axes(fig3);
    b = a;
    fluos = [b.particleFluo3Slice];
    fluos(fluos <= 0) = [];
    fluos(isnan(fluos)) = [];
    fluos = fluos./max(fluos(:));
 if ~isempty(fluos)
    histogram(ax3, fluos, 'Normalization', 'pdf', 'facealpha', .5)
    pd = fitdist(fluos','Lognormal');
    hold on
    x_values = 0:1E-3:1;
    y = pdf(pd,x_values);
    plot(ax3, x_values,y,'-','LineWidth',1)
    hold on
end
xlabel('spot fluo normalized to max')
ylabel('pdf')

coeffText = getDistributionText(pd);
title(["spot fluorescence fit to lognormal";coeffText'])

% paperize(gcf, 4,4)

%%
%%
fig4= figure;
ax4 = axes(fig4);
b = a;
fluos = [];
for k = 1:length(b)
    fluos = [fluos, max(b(k).particleFluo3Slice)];
end  
% fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
fluos = fluos + 100;
%  fluos = fluos./max(fluos(:));
 if ~isempty(fluos)
    histogram(ax4, fluos, 'Normalization', 'pdf', 'facealpha', .5)
    hold on;
    pd = fitdist(fluos','Lognormal');
    x_values = .01:.01:3000;
    y = pdf(pd,x_values);
    plot(ax4, x_values,y,'-','LineWidth',1.5)
 end
 coeffText = getDistributionText(pd);

title(['max spot fluorescence. Lognormal fit';coeffText'])
ylabel('pdf')
xlabel('max fluo + 100 (au)')

% mu = log(mean(fluos)) - 
% varx = var(fluos)  
% sigma2 = log( (1/2) + (1/2)*exp(  -2*mu*sqrt(  exp(4*mu) + 4*varx*exp(2*mu)  )  ) );

%%
fig5= figure;
ax5 = axes(fig5);
b = a;
fluos = [];
for k = 1:length(b)
    fluos = [fluos, mean(b(k).particleFluo3Slice)];
end  
fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
%  fluos = fluos./max(fluos(:));
 if ~isempty(fluos)
    histogram(ax5, fluos, 'Normalization', 'pdf', 'facealpha', .5)
    hold on;
    pd = fitdist(fluos','Lognormal');
    x_values = 1:1:1000;
    y = pdf(pd,x_values);
    plot(ax5, x_values,y,'-','LineWidth',.5)
 end
 coeffText = getDistributionText(pd);

title(['mean spot fluorescence. Lognormal fit';coeffText'])
ylabel('pdf')


%%
figure;
b = a;
fluos = [];
for k = 1:length(b)
    fluos = [fluos, max(b(k).particleFluo3Slice)];
end  
% fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
fluos = log(fluos+100);
%  fluos = fluos./max(fluos(:));
 if ~isempty(fluos)
    h = histogram(fluos, 'Normalization', 'pdf', 'facealpha', .5)
    hold on;
    pd = fitdist(fluos','Normal');
    L = log(54 + 100); %54 aus from Simon's figure. 
    [norm_trunc, phat, phat_ci]  = fitdist_ntrunc(fluos', [L, Inf]);
    x_values = .01:.01:9;
    y = pdf(pd,x_values);
    plot(x_values,y,'-','LineWidth',1)
    plot(x_values,norm_trunc(x_values , phat(1), phat(2)),'-','LineWidth',2)
 end
 coeffText = getDistributionText(pd);

title(['log(max spot fluorescence + 100). Normal fit';coeffText'])
legend('data', 'normal', 'truncated normal')
ylabel('pdf')


%%
figure;
b = a;
fluos = [];
for k = 1:length(b)
    fluos = [fluos, max(b(k).particleFluo3Slice)];
end  
% fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
fluos = log(fluos+100);
[f,xi] = ksdensity(fluos); %estimate kernel density for help finding mode of the distribution
[~, mi] = max(f); %here's the index of the mode
pk = xi(mi); %and here's the fluo value of the mode
pk = 5.6;
fluosWhole = fluos;
fluos(fluos < pk) = [];
%  fluos = fluos./max(fluos(:));
 if ~isempty(fluos)
    yyaxis left;
%     h = histogram(fluos, 'Normalization', 'pdf', 'facealpha', .5);
    hold on;
    pd = fitdist(fluos','HalfNormal', 'mu', pk);
    L = log(54 + 100); %54 aus from Simon's figure. 
%     [norm_trunc, phat, phat_ci]  = fitdist_ntrunc(fluos', [L, Inf]);
    x_values = .01:.01:9;
    y = pdf(pd,x_values);
    plot(x_values,y,'-','LineWidth',2)
    ylim([0, max(y)]);
    yyaxis right;
    h2 = histogram(fluosWhole, 'Normalization', 'pdf', 'facealpha', .5);
    ylim([0, max(h2.Values)]);
    yyaxis left;
    n = makedist('Normal', 'mu', pd.mu, 'sigma', pd.sigma);
    plot(x_values, pdf(n, x_values),'-','LineWidth',2)



%     plot(x_values,norm_trunc(x_values , phat(1), phat(2)),'-','LineWidth',2)
 end
 coeffText = getDistributionText(pd);

title(['log(max spot fluorescence + 100). Normal fit';coeffText'])
legend('data', 'normal', 'truncated normal')
ylabel('pdf')


%%
figure;
b = a;
fluos = [];
pk = 5.5;
for k = 1:length(b)
    fluos = [fluos, max(b(k).particleFluo3Slice)];
end  
% fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
fluos = log(fluos+100);
%  fluos = fluos./max(fluos(:));
 if ~isempty(fluos)
    h = histogram(fluos, 'Normalization', 'pdf', 'facealpha', .5);
    hold on;
    yy = h.Values;
%     x = h.BinEdges(2:end); 
    x = mean([h.BinEdges(1:end-1);h.BinEdges(2:end)]);
    yyhalf = yy(x>pk);
    xhalf = x(x > pk);
%     figure; bar(xhalf, yyhalf);
    yyhalfbigger = [fliplr(yyhalf), yyhalf];
    xhalfbigger = xhalf;
    for k = 1:length(yyhalf)
        xhalfbigger = [xhalf(1)-(k*h.BinWidth),xhalfbigger];
    end
    bar(xhalfbigger, yyhalfbigger, 'facealpha', .5);
    [mu, sigma] = fitnormal(xhalfbigger, yyhalfbigger);
%     yyhalfbigger = [fliplr(yyhalf), yyhalf];
%     pd = fitdist(fluos','Normal');
    pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
    L = log(54 + 100); %54 aus from Simon's figure. 
%     [norm_trunc, phat, phat_ci]  = fitdist_ntrunc(fluos', [L, Inf]);
    x_values = .01:.01:9;
    y = pdf(pd,x_values);
    plot(x_values,y,'-','LineWidth',1)
%     plot(x_values,norm_trunc(x_values , phat(1), phat(2)),'-','LineWidth',2)
    xline(L, 'LineWidth', 2);
 end
 coeffText = getDistributionText(pd);

title(['log(max spot fluorescence + 100). Normal fit';coeffText'])
legend('half data', 'reflected data', 'normal fit to reflected')
ylabel('pdf')


%%
CompiledParticles = CompiledParticles{1};

%%

for p =1:length(CompiledParticles)
    particleFluo = CompiledParticles(p).Fluo;
    particleTime = CompiledParticles(p).Frame;
    plot(particleTime,particleFluo,'ro-')
    waitforbuttonpress
end

%%
% data from S:\Jonathan\Dropbox\ZeldaLizHernan\Data\SingleParticleFitsCompiled
tracesData = fits.other;
rateDataRaw =  fits.rate;
deltaT = fits.t_steady - fits.t_on;
maxFluos = [];
for i = 1:length(tracesData)
    traceFluo = tracesData(i).MS2_raw;
    maxFluos(i) = prctile(smooth(traceFluo),90);%max(smooth(traceFluo,3));
end
rateData = rateDataRaw./deltaT;
figure
hold on
plot(maxFluos,rateData,'o','MarkerFaceColor','r','MarkerEdgeColor','none')
ylabel('RNAP loading rate (a.u)/min')
xlabel('maximum trace fluorescence')

%fit
mdl = fitlm(maxFluos,rateData);
fittedCoeff = mdl.Coefficients;
slopeMean = fittedCoeff{2,1};
slopeError = fittedCoeff{2,2};
XInterceptMean = fittedCoeff{1,1};
XInterceptError = fittedCoeff{1,2};

% show fit on top
plot(maxFluos,XInterceptMean + (maxFluos.*slopeMean),'k-','LineWIdth',2)

% now through the origin
mdl2 = fitlm(maxFluos,rateData,'Intercept',false);
fittedCoeff2 = mdl2.Coefficients;
slopeMean2 = fittedCoeff2{1,1};

% show fit on top
plot(maxFluos,(maxFluos.*slopeMean2),'g-','LineWIdth',2)
hold off
%legend({'data','with intercept, RMSE=136','without intercept, RMSE=156'})
% 
% % a quadratic model
% mdl3 = fitlm((maxFluos).^2,rateData)%,'Intercept',false);
% fittedCoeff3 = mdl3.Coefficients;
% slopeMean3 = fittedCoeff3{1,1};
% figure
% hold on
% plot(sqrt(maxFluos),rateData,'ko')
% plot(sqrt(maxFluos),(sqrt(maxFluos).*slopeMean3),'m-','LineWIdth',2)

%% get average max fluo per AP bin
% add meanAP position and max fluo to struct with particles

% to take errors across embryos, add info about prefixes
positions = fits.MeanAP;
UniqueEmbryos = unique(fits.Prefix);
UniqueEmbryoNames = fits.Prefix;
UniqueEmbryoIDs = [];
for e = 1:length(UniqueEmbryos)
    EmbryoName = UniqueEmbryos{e};
    UniqueEmbryoIDs(strcmpi(UniqueEmbryoNames,EmbryoName)) = e;
end

% add data and info to struct with particles
for i = 1:length(tracesData)
    tracesData(i).t_elon  = fits.t_steady(i)-fits.t_on(i);
    tracesData(i).rate = fits.rate(i);
    tracesData(i).meanAP = positions(i);
    traceFluo = tracesData(i).MS2_raw;
    tracesData(i).maxFluo = prctile(smooth(traceFluo),90);
    tracesData(i).EmbryoID = UniqueEmbryoIDs(i);
    tracesData(i).maxFluoBin = [];
end

%% bin by max fluo
maxFluo = 2700;
Nbins = 30;
maxFluoBins = linspace(0,maxFluo,Nbins);
maxFluoBins = BinData([tracesData.maxFluo],maxFluoBins);
for i = 1:length(tracesData)
    tracesData(i).maxFluoBin = maxFluoBins(i);
end

% calculate average rate per max fluo bin
allData = nan(Nbins,50);
meanRatePerMaxFluoBin = nan(1,Nbins);
for i = 1:Nbins
    tracesThisBin = tracesData([tracesData.maxFluoBin] == i);
    meanRatePerMaxFluoBin(i) = mean([tracesThisBin.rate]);
    SEMRatePerMaxFluoBin(i) = std([tracesThisBin.rate])./length([tracesThisBin.rate])-1;
    allData(i,1:length([tracesThisBin.rate])) = [tracesThisBin.rate];
end

Nembryos = sum(~isnan(allData),2);
figure(1)
hold on
%plot(1:Nbins,meanRatePerMaxFluoBin,'o')
%errorbar(1:Nbins,meanRatePerMaxFluoBin,SEMRatePerMaxFluoBin,'r','CapSize',0)
plot(1:Nbins,allData,'o','MarkerEdgeColor','none','MarkerFaceColor','k')
errorbar(1:Nbins,nanmean(allData,2),nanstd(allData,[],2)./(Nembryos-1),'ro','CapSize',0,...
    'MarkerFaceColor','r','MarkerEdgeColor','none','LineWidth',2)
%fit with intercept
mdl = fitlm(1:Nbins,meanRatePerMaxFluoBin);
fittedCoeff = mdl.Coefficients;
slopeMean = fittedCoeff{2,1};
Yint = fittedCoeff{1,1};
plot(1:Nbins,Yint+([1:Nbins].*slopeMean),'g-')
%fit without intercept
mdl2 = fitlm(1:Nbins,meanRatePerMaxFluoBin,'Intercept',false);
fittedCoeff2 = mdl2.Coefficients;
slopeMean2 = fittedCoeff2{1,1};
plot(1:Nbins,[1:Nbins].*slopeMean2,'m-')
xlabel('max fluorescence bin')
ylabel('mean rate across spots in the bin')
legend(['R^2 with intercept= ' num2str(mdl.Rsquared.Ordinary)],...
    ['R^2 without intercept= ' num2str(mdl2.Rsquared.Ordinary)])
hold off

figure(2)
plot(1:Nbins,allData,'o','MarkerEdgeColor','none','MarkerFaceColor','k')
errorbar(1:Nbins,nanmean(allData,2),nanstd(allData,[],2)./(Nembryos-1),'ro','CapSize',0,...
    'MarkerFaceColor','r','MarkerEdgeColor','none','LineWidth',2)
mdl = fitlm(1:Nbins,meanRatePerMaxFluoBin);
fittedCoeff = mdl.Coefficients;
slopeMean = fittedCoeff{2,1};
Yint = fittedCoeff{1,1};
plot(1:Nbins,Yint+([1:Nbins].*slopeMean),'k-')
xlabel('max fluorescence bin')
ylabel('mean rate across spots in the bin')
hold off



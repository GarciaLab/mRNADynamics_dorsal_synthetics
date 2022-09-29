function P2PCalibration

Prefixes = {'2020-10-14-P2Pms2v1-2xDl_syntheticssettings_1','2020-10-14-P2Pms2v1-2xDl_syntheticssettings_2',...
    '2020-10-14-P2Pms2v1-2xDl_syntheticssettings_3','2020-10-14-P2Pms2v1-2xDl_syntheticssettings_4',...
    '2020-10-14-P2Pms2v1-2xDl_syntheticssettings_5','2020-10-14-P2Pms2v1-2xDl_syntheticssettings_6'};
% for p =1:length(Prefixes)
%     Prefix = Prefixes{p}; 
%     AddParticlePosition(Prefix,'SkipAll','ApproveAll')
% end
% Sum particle fluorescence per AP bin

ResultsFolder = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox';
APbinIDs = linspace(0,1,21);
AccFluoPerBin = zeros(length(Prefixes),length(APbinIDs));
for prfx = 1:length(Prefixes)

    PrefixAccumulatedFluoPerBin = zeros(size(APbinIDs));
    resultsPath = [ResultsFolder '\' Prefixes{prfx}];
    load([resultsPath '\CompiledParticles.mat'],'CompiledParticles')
    CompiledParticles = CompiledParticles{1};
    load([resultsPath '\FrameInfo.mat']);
    PrefixFrames = 1:length(FrameInfo);
    PrefixAbsTime = [FrameInfo.Time]./60;
    
    NCs = [CompiledParticles.nc];
    particleIdxs = 1:length(CompiledParticles);
    NC13ParticleIdx =particleIdxs(NCs==13);
    blankParticleFluo = zeros(size(PrefixFrames));
    
    for prt = NC13ParticleIdx
        particleFrames = CompiledParticles(prt).Frame;
        particleFluo = CompiledParticles(prt).Fluo3;
        blankParticleFluo(particleFrames) = particleFluo;
        %plot(PrefixAbsTime,blankParticleFluo,'ro-')
        %waitforbuttonpress
        particleTotalmRNA1 = trapz(PrefixAbsTime,blankParticleFluo);
        
        particleAPpos = CompiledParticles(prt).MeanAP;
        particleTotalmRNA2 = CompiledParticles(prt).TotalmRNA;
        
        %[particleTotalmRNA1 particleTotalmRNA2]
        
        [~,APbin] = min(abs(APbinIDs-particleAPpos));
        if ~isempty(particleTotalmRNA1)
            PrefixAccumulatedFluoPerBin(APbin) = PrefixAccumulatedFluoPerBin(APbin)+particleTotalmRNA1;
        end
    end
    AccFluoPerBin(prfx,:) = PrefixAccumulatedFluoPerBin;
end

% sum the number of nuclei per AP bin

%these are the frames in which we count the number of nuclei
StableFrames = [56,54,45,84,28,60];
ResultsFolder = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox';
APbinIDs = linspace(0,1,21);
NumberOfNucleiPerBin = zeros(length(Prefixes),length(APbinIDs));
for prfx = 1:length(Prefixes)
    PrefixNumberOfNucleiPerBin = zeros(size(APbinIDs));
    resultsPath = [ResultsFolder '\' Prefixes{prfx}];
    load([resultsPath '\' Prefixes{prfx} '_lin.mat'])
    
    NCs = [schnitzcells.cycle];
    NucIdxs = 1:length(schnitzcells);
    NC13NucIdx =NucIdxs(NCs==13);
    
    for nuc = NC13NucIdx
        nucFrames = schnitzcells(nuc).frames;
        if ismember(StableFrames(prfx),nucFrames)
            nucAPpos = schnitzcells(nuc).APPos;
            nucAPpos = nucAPpos(nucFrames==StableFrames(prfx));
            [~,APbin] = min(abs(APbinIDs-nucAPpos));
            if ~isempty(nucAPpos)
                PrefixNumberOfNucleiPerBin(APbin) = PrefixNumberOfNucleiPerBin(APbin) + 1;
            end
        end
    end
    NumberOfNucleiPerBin(prfx,:)=PrefixNumberOfNucleiPerBin;
end


% Plot MS2 results
figure
minNucPerBin = 3;
NumberOfNucleiPerBin(NumberOfNucleiPerBin<minNucPerBin) = nan;
NumberOfNucleiPerBin(NumberOfNucleiPerBin==0) = nan;

hold on
YDataMS2 = AccFluoPerBin./NumberOfNucleiPerBin;
meanY = nanmean(YDataMS2);
NY = sum(~isnan(YDataMS2));
seY = nanstd(YDataMS2)./NY;
errorbar(APbinIDs,meanY,seY,'k','CapSize',0)
plot(APbinIDs,YDataMS2,'ko')
ylabel('accumulated mRNA (AU \times frame)')
xlabel('AP position (x/L)')

% Load FISH data
FISHPath = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox\mRNA_FISH_calibration_data/FISH mRNA per Nucleus vs AP Position.csv'

FISHData = readcell(FISHPath);

APPosFISH = cell2mat(FISHData((2:end),1));
MeanFISH = cell2mat(FISHData((2:end),2));

figure
plot(APPosFISH,MeanFISH,'-o','Color',[0.1 0.8 0.2],'LineWIdth',2)
ylabel('number of mRNA per nc13 nucleus')
xlabel('AP position (x/L)')
legend('smFISH data from Garcia 2013')


%% plot side by side
%close all
clearvars -except APbinIDs YDataMS2 APPosFISH MeanFISH NY
meanYDataMS2 = nanmean(YDataMS2);

meanYDataMS2 = meanYDataMS2(7:end); % more anterior bins are sketch right now
APbinIDs2 = APbinIDs(7:end);
%NY = NY(7:end);
seY = nanstd(YDataMS2)./NY;
seY = seY(7:end);

figure
yyaxis left
%plot(APbinIDs2,meanYDataMS2,'o','MarkerFaceColor','w','MarkerEdgeColor','k','LineWIdth',1.5)
errorbar(APbinIDs2,meanYDataMS2,seY,'o','MarkerFaceColor','w','MarkerEdgeColor','k','CapSize',0)
ylabel('accumulated mRNA (AU)')

yyaxis right
plot(APPosFISH,MeanFISH,'o','MarkerEdgeColor',[0.1 0.8 0.2],'MarkerFaceColor','w')%,'LineWIdth',1.5)
ylabel('number of mRNAs (smFISH)')


downsampledFISH_Y = [];
APbinIDs_downsampledFISH = [];
finalMS2_Y = [];
APbinIDs_finalMS2 = [];
for b = 1:length(APbinIDs2)
    APValue = APbinIDs2(b);
    [diff,FISHidx] = min(abs((APPosFISH-APValue)));
    if diff < 0.05
        MS2Value = meanYDataMS2(b);
        FISHValue = MeanFISH(FISHidx);       
        finalMS2_Y = [finalMS2_Y MS2Value];
        APbinIDs_finalMS2 = [APbinIDs_finalMS2 APValue];
        downsampledFISH_Y = [downsampledFISH_Y FISHValue];
        APbinIDs_downsampledFISH = [APbinIDs_downsampledFISH APValue];
    end
end

figure
yyaxis right
plot(APbinIDs_finalMS2,finalMS2_Y,'o','MarkerFaceColor','w','MarkerEdgeColor','k','LineWIdth',1.5)
ylabel('accumulated mRNA (AU)')

yyaxis left
plot(APbinIDs_downsampledFISH,downsampledFISH_Y,'*','MarkerEdgeColor',[0.1 0.8 0.2],'MarkerFaceColor','w','LineWIdth',1.5)
ylabel ('downsampled FISH')

%% Calibrate with y intercept

% now do the fitting allowing a y intercept
XdataVector = [finalMS2_Y]';%
YdataVector = [downsampledFISH_Y]';
%FluoError = [SDVector;SDVector];

regX = [ones(size(XdataVector)) XdataVector];% the padding with ones is needed for the thing to run
[b,bint,r,rint,stats] = regress(YdataVector,regX);
slope1 = b(2);
minslope = bint(2,2) ; % top of 95% confidence interval
maxslope = bint (2,1);% bottom of 95% confidence interval
Yintersect = b(1);
Rsquared = stats(1);
% add a zero to the X data for plotting purposes
XdataVector = [XdataVector];

figure
hold on
plot(XdataVector,Yintersect+XdataVector*slope1,'k-','LineWidth',1)
plot(finalMS2_Y,downsampledFISH_Y,'ko')
xlabel('produced mRNA (MS2)')
ylabel('number of mRNA (smFISH)')
title(['R^2 = ' num2str(Rsquared) ' ; slope = ' num2str(slope1) '\pm ' num2str(abs(slope1-minslope))])
hold off


%% calibrate through origin

% now do the fitting forcing it to go through zero
XdataVector = finalMS2_Y;%
YdataVector = downsampledFISH_Y;

mdl = fitlm(XdataVector,YdataVector,'Intercept',false);
slope2 = mdl.Coefficients.Estimate;
slopeError = mdl.Coefficients.SE;
%yfit = linspace(0,750,20).*slope;
Rsqrd = mdl.Rsquared.Ordinary;

yfit = XdataVector*slope2; % calculate fitted line

% % calculate R2          
% f = yfit;
% Bbar = mean(XdataVector');
% SStot = sum((XdataVector' - Bbar).^2);
% SSreg = sum((f - Bbar).^2);
% SSres = sum((XdataVector' - f).^2);
% R2 = 1 - SSres/SStot;
% R = corrcoef(YdataVector,XdataVector');
% Rsq = R(1,2).^2;

% add a zero to data for plotting purposes
XdataVector = [0 XdataVector];
yfit = [0 yfit];

figure
hold on
plot(finalMS2_Y,downsampledFISH_Y,'ko')
plot(XdataVector,yfit,'k-','LineWidth',1)
hold off
title(['Slope = ' num2str(slope2) '\pm' num2str(slopeError) ' ; R^2 = ' num2str(Rsqrd)])
xlabel('produced mRNA (MS2)')
ylabel('number of mRNA (smFISH)')

%% Calculate the fluorescence of 1RNAP

a = 1/slope2;
alphaE = slopeError;
Velo = 1.5; %elongation speed in Kbp/min
VeloE = 0.14;
L = 5.2; %gene length in Kbp
FluoRNAP = (Velo*a)*(1/L); %AUs per RNAP

% now propagate the error in the constants
X1 = a;
E1 = VeloE;
X2 = Velo;
E2 = alphaE;

% propagate the error
FluoRNAPError = ((X1*X2) * sqrt((E1/X1)^2 + (E2/X2)^2))*(1/L);

percentError = FluoRNAPError/FluoRNAP;

%% Compare results with spot intensities
clearvars -except FluoRNAP FluoRNAPError
close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])%, 'dorsalResultsDatabase')
a = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);


%% all spots combined, the minimum fluorescence per trace
b = a;
fluos = [];
fluoError = [];
NucIdx = [];
for k = 1:length(b)
    thisFluoTrace = b(k).particleFluo3Slice;
    thisTraceDimmest = nanmean(thisFluoTrace(thisFluoTrace>0));
%     thisTraceSorted = sort(thisFluoTrace);
%     thisTraceDimmest = thisTraceSorted(1:min(4,length(thisTraceSorted))); 
%     thisTraceDimmest(thisTraceDimmest<0) = [];
    if ~isempty(thisTraceDimmest)
        if ~isnan(thisTraceDimmest)
            fluos = [fluos mean(thisTraceDimmest)];
            fluoError = [fluoError b(k).particleFluoError];
            NucIdx = [NucIdx k];
        end
    end
end  

fluos3 = log(fluos);
absFluos = log(fluos./FluoRNAP);
% multiply error by sqrt(3) because we're using 3 slice fluo values
absFluoError3 = log((fluoError./FluoRNAP)*sqrt(3));
fluoError3 = log((fluoError*sqrt(3)));

figure
histogram(fluos3, 'Normalization', 'probability', 'facealpha', .5)
hold on
histogram(fluoError3, 'Normalization', 'probability', 'facealpha', .5)
legend('weak spots','fluo error x3')
title('min spot fluorescence. Lognormal fit')
ylabel('probability')
xlabel('log(fluo)(AU)')
xlim([1 7])
hold off

figure
% plot the histograms in RNAP scale
histogram(absFluos, 'Normalization', 'probability', 'facealpha', .5)
hold on
histogram(absFluoError3, 'Normalization', 'probability', 'facealpha', .5)
ylabel('probability')
xlabel('log(RNAP)')
legend('weak spots','fluo error x3')
hold off
xlim(log([exp(1)/FluoRNAP exp(7)/FluoRNAP]))



figure
scatter(log(fluos),log(fluoError))
xlabel('spot fluo')
ylabel('fluo error')


%%

allTogether = [fluos;fluoError;NucIdx]';
allTogetherSorted = sortrows(allTogether,2);

nData = 455;
dimmestNucleiIdx = allTogetherSorted(1:nData,3);
fluos2 = [];
fluos3=[];
errors = [];

for n = 1:length(dimmestNucleiIdx)
    idx = dimmestNucleiIdx(n)
    trace = b(idx).particleFluo3Slice;
    traceDimmestFrames = sort(trace);
    traceDimmestFrames = traceDimmestFrames(1:min(3,length(trace)));
    traceError = b(idx).particleFluoError;
    
    fluos2 = [fluos2 traceDimmestFrames];
    fluos3 = [fluos3 nanmean(traceDimmestFrames)];
    errors = [errors traceError];
end

figure
hold on
histogram(fluos2,'Normalization', 'pdf', 'facealpha', .5)
histogram(errors(errors>0)*sqrt(3),'Normalization', 'pdf', 'facealpha', .5)
hold off
legend('fluo','error')

figure
scatter(fluos3,errors*sqrt(3))
xlabel('spot fluo')
ylabel('error')

% %
% fig3 = figure;
% ax3 = axes(fig3);
%     b = a;
%     fluos = [b.particleFluo3Slice];
%     fluos(fluos <= 0) = [];
%     fluos(isnan(fluos)) = [];
%     fluos = fluos./max(fluos(:));
%  if ~isempty(fluos)
%     histogram(ax3, fluos, 'Normalization', 'pdf', 'facealpha', .5)
%     pd = fitdist(fluos','Lognormal');
%     hold on
%     x_values = 0:1E-3:1;
%     y = pdf(pd,x_values);
%     plot(ax3, x_values,y,'-','LineWidth',1)
%     hold on
% end
% xlabel('spot fluo normalized to max')
% ylabel('pdf')
% coeffText = getDistributionText(pd);
% title(["spot fluorescence fit to lognormal";coeffText'])
% paperize(gcf, 4,4)




end






    
    





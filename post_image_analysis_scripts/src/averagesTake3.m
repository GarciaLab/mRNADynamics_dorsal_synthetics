function [embryoRNA embryoRNAError embryoMetric embryoMetricError] = ...
    averagesTake3(DataType,NotEnhancerName,numBins,metric,fiducialTime,errorgroup,Color,ax)
% USES ABSOLUTE DV POSITION

% metric can be 'maxfluo', 'accumulatedfluo' or 'fraction'
% errorgroup is over what the error is taken, 'embryos' or 'nuclei'. For
% fraction the error over nuclei is bootstraped.

% load everything
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

% add absolute DV
combinedCompiledProjects_allEnhancers = mapFluoToDV(combinedCompiledProjects_allEnhancers);

% we'll use the struct called 'combinedCompiledProjects_allEnhancers', wich
% contains one entry per nucleus.
% First, make a smaller struct called 'enhancerStruct'containing only the nc12 nuclei for the enhancer
% specified by the 'DataTye' argument.
for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

% define rules for what data we'll consider in the analysis and plots
dataTypeIncludeRule = contains({combinedCompiledProjects_allEnhancers.dataSet},DataType);
dataTypeExcludeRule = ones(1,length(combinedCompiledProjects_allEnhancers));
for i = 1:length(NotEnhancerName)
    excludeString = NotEnhancerName{i};
    Exclude = ~contains(lower({combinedCompiledProjects_allEnhancers.dataSet}),lower(excludeString));
    dataTypeExcludeRule = dataTypeExcludeRule & Exclude;
end
dataTypeRule = dataTypeIncludeRule & dataTypeExcludeRule;

% tempRule = ({combinedCompiledProjects_allEnhancers.dataSet} == "1Dg11_2xDl" |...
%       {combinedCompiledProjects_allEnhancers.dataSet} == "1Dg11_FFF"  );
% dataTypeRule = dataTypeRule & tempRule;

enhancerStruct = combinedCompiledProjects_allEnhancers(dataTypeRule &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

% bin nuclei 

% ****IMPORTANT!!!**** this line uses the standard dorsal fluorescence, not the
% arbitrary one at a given 'fiducial time'
nucleiPositions = [enhancerStruct.AbsDV];
binValues = linspace(0,0.66,numBins);
binnedNuclearDVPos = BinData(nucleiPositions,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).DVpos = binnedNuclearDVPos(n);
end

coveredBins = unique([enhancerStruct.DVpos]);
binValues = binValues(coveredBins);

%% bin nuclei in the same way we do in getEmbryoDataForFits
dorsalVals = linspace(0,0.66,numBins);
% binSize = diff(dorsalVals(1:2));
% for n = 1:length(enhancerStruct)
%     FluoFeature = enhancerStruct(n).dorsalFluoFeature;
%     bin = find(FluoFeature-dorsalVals<0,1,'first')-1;
%     enhancerStruct(n).dorsalFluoBin3 = bin;
% end
% coveredBins = unique([enhancerStruct.dorsalFluoBin3]);
% binValues = dorsalVals(1:end-1) + (binSize/2);

%% now make one struct per embryo per bin 
% define some filters
minEmbryosPerBin = 3;
minNucleiPerEmbryoPerBin = 2;
minOnset = 2; % (min) earliest possible spot detection time to be counted
maxOnset = 8; %(min) latest possible spot detection time to be counted
maxMaxFluo = 1500; % (a.u)

%store everything in these arrays
mean_maxFluo_acrossNuclei_perBin = [];
se_maxFluo_acrossNuclei_perBin = [];
mean_accFluo_acrossNuclei_perBin = [];
se_accFluo_acrossNuclei_perBin = [];
mean_fraction_acrossNuclei_perBin = [];
%
mean_maxFluo_acrossEmbryos_perBin = [];
se_maxFluo_acrossEmbryos_perBin = [];
mean_accFluo_acrossEmbryos_perBin = [];
se_accFluo_acrossEmbryos_perBin = [];
mean_fraction_acrossEmbryos_perBin = [];
se_fraction_acrossEmbryos_perBin = [];

numEmbryos = 60; %some arbitrarily large number
fractionPerEmbryoPerBin = nan(numEmbryos,length(coveredBins));
fractionPerEmbryoPerBin_cellarray = {};

OnNucleiPerBin = [];
OffNucleiPerBin = [];
for b = 1:length(coveredBins)
    
    %binID = coveredBins(b);
    %binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin2]== binID);
    binStruct = enhancerStruct([enhancerStruct.DVpos]== b);

    if ~isempty(binStruct)
    
        activeNuc_Bin = length([binStruct.particleTimeOn]);
        inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;
        % take the means and errors across nuclei first because it's easy 
        %filter spots that are way too bright to be real spots
        particlesMaxFluos = [binStruct.particleFluo95];
        particlesMaxFluos = particlesMaxFluos(particlesMaxFluos<maxMaxFluo);
        particlesMaxFluos = particlesMaxFluos(particlesMaxFluos>0);
        mean_maxFluo_acrossNuclei_perBin(b) = nanmean(particlesMaxFluos);
        se_maxFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleFluo95])./sqrt(activeNuc_Bin);
        mean_accFluo_acrossNuclei_perBin(b) =  nanmean([binStruct.particleAccumulatedFluo]);
        se_accFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleAccumulatedFluo])./sqrt(activeNuc_Bin);
        mean_fraction_acrossNuclei_perBin(b) = activeNuc_Bin/length(binStruct);
        %filter spurious time ons due to errors
        particlesTimeOns = [binStruct.particleTimeOn];
        particlesTimeOns = particlesTimeOns(particlesTimeOns>minOnset);
        particlesTimeOns = particlesTimeOns(particlesTimeOns<maxOnset);
        mean_timeOn_acrossNuclei_perBin(b) = nanmean(particlesTimeOns);
        se_timeOn_acrossNuclei_perBin(b) = nanstd(particlesTimeOns)./sqrt(activeNuc_Bin);
        mean_duration_acrossNuclei_perBin(b) = nanmean([binStruct.particleDuration]);
        se_duration_acrossNuclei_perBin(b) = nanstd([binStruct.particleDuration])./sqrt(activeNuc_Bin);    
        totalRNA_acrossNuclei_perBin(b) = nansum([binStruct.particleAccumulatedFluo]);

        % save this for bootstrapping later
        OnNucleiPerBin(b) = activeNuc_Bin;
        OffNucleiPerBin(b) = inactiveNuc_Bin;

        % now deal with the mean and errors across embryos
        [uniqueEmbryos, ~, J]=unique({binStruct.prefix});
        occurences = histc(J, 1:numel(uniqueEmbryos));
        numEmbryos = length(occurences);
        qualityEmbryos = occurences >= minNucleiPerEmbryoPerBin;
        maxFluo_perEmbryo = [];
        accFluo_perEmbryo = [];
        fraction_perEmbryo = []; 
        nuclei_perEmbryo = [];
        timeOn_perEmbryo = [];
        duration_perEmbryo = [];
        totalRNA_perEmbryo = [];

        if numEmbryos >= minEmbryosPerBin     
            for e = 1:numEmbryos
                embryoPrefix = uniqueEmbryos{e};
    %            if qualityEmbryos(e) %if this prefix has enough nuclei
                embryoStruct = binStruct(strcmpi({binStruct.prefix},embryoPrefix));
                nuclei_perEmbryo(e) = length(embryoStruct);
                % filter out spots that are too bright
                particlesMaxFluos = [embryoStruct.particleFluo95];
                particlesMaxFluos = particlesMaxFluos(particlesMaxFluos<maxMaxFluo);
                particlesMaxFluos = particlesMaxFluos(particlesMaxFluos>0);
                maxFluo_perEmbryo(e) = nanmean(particlesMaxFluos);
                accFluo_perEmbryo(e) =  nanmean([embryoStruct.particleAccumulatedFluo]);
                fraction_perEmbryo(e) = length([embryoStruct.particleTimeOn])/length(embryoStruct);
                % clean up the time ons to filter out outliers from errors
                perEmbryoTimeOns = [embryoStruct.particleTimeOn];
                perEmbryoTimeOns = perEmbryoTimeOns(perEmbryoTimeOns>minOnset);
                perEmbryoTimeOns = perEmbryoTimeOns(perEmbryoTimeOns<maxOnset);
                timeOn_perEmbryo(e) = nanmean(perEmbryoTimeOns);
                duration_perEmbryo(e) = nanmean([embryoStruct.particleDuration]);
                totalRNA_perEmbryo(e) =  nansum([embryoStruct.particleAccumulatedFluo]).*fraction_perEmbryo(e);
    %            end
            end
        end
    
        %add things taken across embryos to the arrays where we store stuff
        mean_maxFluo_acrossEmbryos_perBin(b) = nanmean(maxFluo_perEmbryo);
        se_maxFluo_acrossEmbryos_perBin(b) = nanstd(maxFluo_perEmbryo)./sqrt(numEmbryos);
        mean_accFluo_acrossEmbryos_perBin(b) = nanmean(accFluo_perEmbryo);
        se_accFluo_acrossEmbryos_perBin(b) = nanstd(accFluo_perEmbryo)./sqrt(numEmbryos);
        mean_fraction_acrossEmbryos_perBin(b) = nanmean(fraction_perEmbryo);
        se_fraction_acrossEmbryos_perBin(b) = std(fraction_perEmbryo)./sqrt(numEmbryos);
        mean_timeOn_acrossEmbryos_perBin(b) = nanmean(timeOn_perEmbryo);
        se_timeOn_acrossEmbryos_perBin(b) = nanstd(timeOn_perEmbryo)./sqrt(numEmbryos);
        mean_duration_acrossEmbryos_perBin(b) = nanmean(duration_perEmbryo);
        se_duration_acrossEmbryos_perBin(b) = nanstd(duration_perEmbryo)./sqrt(numEmbryos);
        mean_totalRNA_acrossEmbryos_perBin(b) = nanmean(totalRNA_perEmbryo);
        se_totalRNA_acrossEmbryos_perBin(b) = nanstd(totalRNA_perEmbryo)./sqrt(numEmbryos);

        %clear  nuclei_per_Embryo maxFluo_perEmbryo accFluo_perEmbryo fraction_perEmbryo timeOn_perEmbryo

        fractionPerEmbryoPerBin(1:length(fraction_perEmbryo),b) = fraction_perEmbryo;
        fractionPerEmbryoPerBin_cellarray{b} = fraction_perEmbryo;
        
    else % if this enhancer had no embryos with nuclei in this bin
        
        OnNucleiPerBin(b) = nan;
        OffNucleiPerBin(b) = nan;
        
        mean_maxFluo_acrossEmbryos_perBin(b) = nan;
        se_maxFluo_acrossEmbryos_perBin(b) = nan;
        mean_accFluo_acrossEmbryos_perBin(b) = nan;
        se_accFluo_acrossEmbryos_perBin(b) = nan;
        mean_fraction_acrossEmbryos_perBin(b) = nan;
        se_fraction_acrossEmbryos_perBin(b) = nan;
        mean_timeOn_acrossEmbryos_perBin(b) = nan;
        se_timeOn_acrossEmbryos_perBin(b) = nan;
        mean_duration_acrossEmbryos_perBin(b) = nan;
        se_duration_acrossEmbryos_perBin(b) = nan;
        mean_totalRNA_acrossEmbryos_perBin(b) = nan;
        se_totalRNA_acrossEmbryos_perBin(b) = nan;

        
    end
     
end

%% Integrate the mean total mRNA produced per nucleus (inactive and active ones) across dorsal fluo bins
% to show the total mRNA produced per embryo

embryoRNA = nansum(mean_totalRNA_acrossEmbryos_perBin);
embryoRNAError = sqrt(nansum(se_totalRNA_acrossEmbryos_perBin).^2); % propagated error

%% Average the metric across Dorsal fluo bins




%% make figures

if strcmpi(metric,'maxfluo') 
    if strcmpi(errorgroup,'nuclei')
    errorbar(ax,binValues,mean_maxFluo_acrossNuclei_perBin,se_maxFluo_acrossNuclei_perBin,'ro-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
    errorbar(ax,binValues,mean_maxFluo_acrossEmbryos_perBin,se_maxFluo_acrossEmbryos_perBin,'ko-','CapSize',0,'LineWidth',1.5,...
       'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('maximum spot fluorescence')
%set(gca,'XScale','log')
% ylim([0 600])
ylim([0,600])
xlim([0 0.66])
%legend('across nuclei','across embryos')
embryoMetric = nanmean(mean_maxFluo_acrossEmbryos_perBin);
embryoMetricError = nanstd(mean_maxFluo_acrossEmbryos_perBin)./sqrt(sum(~isnan(mean_maxFluo_acrossEmbryos_perBin)));


elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_accFluo_acrossNuclei_perBin,se_accFluo_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
     'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
     elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_accFluo_acrossEmbryos_perBin,se_accFluo_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('accumulated fluorescence')
%set(gca,'XScale','log')
xlim([0 0.66])
ylim([0 1200])
embryoMetric = nanmean(mean_accFluo_acrossEmbryos_perBin);
embryoMetricError = nanstd(mean_accFluo_acrossEmbryos_perBin)./sqrt(sum(~isnan(mean_accFluo_acrossEmbryos_perBin)));



elseif contains(lower(metric),'fraction')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_fraction_acrossNuclei_perBin,btsrp_error_fraction_perBin,'o-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_fraction_acrossEmbryos_perBin,se_fraction_acrossEmbryos_perBin,'ko-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
        
        %plot(ax,binValues,fractionPerEmbryoPerBin,'o','MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    
        %plotSpread(fractionPerEmbryoPerBin_cellarray)
    end
xlabel('Dorsal concentration (AU)')
ylabel('fraction of active nuclei')
%set(gca,'XScale','log')
ylim([0 1.1])
xlim([0 0.66])
embryoMetric = nanmean(mean_fraction_acrossEmbryos_perBin);
embryoMetricError = nanstd(mean_fraction_acrossEmbryos_perBin)./sqrt(sum(~isnan(mean_fraction_acrossEmbryos_perBin)));



elseif contains(lower(metric),'timeon')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_timeOn_acrossNuclei_perBin,se_timeOn_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_timeOn_acrossEmbryos_perBin,se_timeOn_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('turn on time')
%set(gca,'XScale','log')
ylim([0 10])
xlim([0 0.66])
embryoMetric = nanmean(mean_timeOn_acrossEmbryos_perBin);
embryoMetricError = nanstd(mean_timeOn_acrossEmbryos_perBin)./sqrt(sum(~isnan(mean_timeOn_acrossEmbryos_perBin)));



elseif contains(lower(metric),'total')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,totalRNA_acrossNuclei_perBin,se_timeOn_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_totalRNA_acrossEmbryos_perBin,se_totalRNA_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('total produced mRNA')
%ylim([0 1800])
xlim([0 0.66])
%set(gca,'YScale','log')
embryoMetric = nanmean(mean_totalRNA_acrossEmbryos_perBin);
embryoMetricError = nanstd(mean_totalRNA_acrossEmbryos_perBin)./sqrt(sum(~isnan(mean_totalRNA_acrossEmbryos_perBin)));



elseif contains(lower(metric),'duration')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_duration_acrossNuclei_perBin,se_duration_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_duration_acrossEmbryos_perBin,se_duration_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('spot duration (min)')
ylim([0 6])
xlim([0 0.66])
%set(gca,'YScale','log')
embryoMetric = nanmean(mean_duration_acrossEmbryos_perBin);
embryoMetricError = nanmean(mean_duration_acrossEmbryos_perBin)./sqrt(sum(~isnan(mean_duration_acrossEmbryos_perBin)));
    

end
%set(gca,'Fontsize',20);



end 

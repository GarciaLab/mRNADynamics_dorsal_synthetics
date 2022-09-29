function [FractionsPerEmbryo,TimeOnsPerEmbryo, MaxFluoPerEmbryo] =  getPerEmbryoDataForFits(includeString,excludeString,dorsalVals)
% DataTypePattern shluld be a string corresponding to a tab in
% dataStatus.xls

% this function is based on plotgreendata_slim_SA. Check that code first of
% you want to add functionalities such as getting the per embryo values of
% other metrics or getting per nucleus metrics.

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

fiducialTime = 6; %mins

for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

% define rules for what data we'll consider in the analysis and plots
dataTypeIncludeRule = contains({combinedCompiledProjects_allEnhancers.dataSet},includeString);
dataTypeExcludeRule = ~contains(lower({combinedCompiledProjects_allEnhancers.dataSet}),lower(excludeString));
dataTypeRule = dataTypeIncludeRule & dataTypeExcludeRule;
enhancerStruct = combinedCompiledProjects_allEnhancers(dataTypeRule &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));


% bin nuclei

for n = 1:length(enhancerStruct)
    FluoFeature = enhancerStruct(n).dorsalFluoFeature;
    bin = find(FluoFeature-dorsalVals<0,1,'first')-1;
    enhancerStruct(n).dorsalFluoBin3 = bin;
end

% these bins are represented in the data
coveredBins = unique([enhancerStruct.dorsalFluoBin3]);
%binValues = binLimits(coveredBins);

% define some filters
minEmbryosPerBin = 3;
minNucleiPerEmbryoPerBin = 2;
minOnset = 2; % earliest possible spot detection time to be counted
maxOnset = 8; % latest possible spot detection time to be counted
minMaxFluo = 0; % dimmest considered spots
maxMaxFluo = 1500; %brightest considered spots

FractionsPerEmbryo = nan(50,length(dorsalVals)-1);
TimeOnsPerEmbryo = nan(50,length(dorsalVals)-1);
MaxFluoPerEmbryo = nan(50,length(dorsalVals)-1);

% now make one struct per bin
for b = 1:length(dorsalVals)-1
    
    binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin3]== b);
    
    if ~isempty(binStruct)
        activeNuc_Bin = length([binStruct.particleTimeOn]);
        inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;

        % now deal with the mean and errors across embryos
        [uniqueEmbryos, ~, J]=unique({binStruct.prefix});
        occurences = histc(J, 1:numel(uniqueEmbryos));
        numEmbryos = length(occurences);
        qualityEmbryos = occurences >= minNucleiPerEmbryoPerBin;
        fraction_perEmbryo = []; 
        nuclei_perEmbryo = [];
        timeOn_perEmbryo = [];
        maxfluo_perEmbryo = [];

        % and now make one struct per embryo per bin
        if numEmbryos >= minEmbryosPerBin     
            for e = 1:numEmbryos
                embryoPrefix = uniqueEmbryos{e};
                embryoStruct = binStruct(strcmpi({binStruct.prefix},embryoPrefix));
                nuclei_perEmbryo(e) = length(embryoStruct);
                fraction_perEmbryo(e) = length([embryoStruct.particleTimeOn])/length(embryoStruct);
                % filter the onset times to 2<true<8 (min since anaphase)
                perEmbryoTimeOns = [embryoStruct.particleTimeOn];
                perEmbryoTimeOns = perEmbryoTimeOns(perEmbryoTimeOns>minOnset);
                perEmbryoTimeOns = perEmbryoTimeOns(perEmbryoTimeOns<maxOnset);
                timeOn_perEmbryo(e) = nanmean(perEmbryoTimeOns);
                %now do the max fluos and filter too
                perEmbryoMaxFluos = [embryoStruct.particleFluo95];
                perEmbryoMaxFluos = perEmbryoMaxFluos(perEmbryoMaxFluos>minMaxFluo);
                perEmbryoMaxFluos = perEmbryoMaxFluos(perEmbryoMaxFluos<maxMaxFluo);
                maxfluo_perEmbryo(e) = nanmean(perEmbryoMaxFluos);
            end
        end
    else % there weren't embryos in this bin
        fraction_perEmbryo = nan;
        timeOn_perEmbryo = nan;
        maxfluo_perEmbryo = nan;
    end
    
%     mean_fraction_acrossEmbryos_perBin(b) = nanmean(fraction_perEmbryo);
%     se_fraction_acrossEmbryos_perBin(b) = std(fraction_perEmbryo)./sqrt(numEmbryos);
%     mean_timeOn_acrossEmbryos_perBin(b) = nanmean(timeOn_perEmbryo);
%     se_timeOn_acrossEmbryos_perBin(b) = nanstd(timeOn_perEmbryo)./sqrt(numEmbryos);

    %save the single embryo datum for later use
    FractionsPerEmbryo(1:length(fraction_perEmbryo),b) = fraction_perEmbryo;
    TimeOnsPerEmbryo(1:length(timeOn_perEmbryo),b) = timeOn_perEmbryo;
    MaxFluoPerEmbryo(1:length(maxfluo_perEmbryo),b) = maxfluo_perEmbryo;
     
end

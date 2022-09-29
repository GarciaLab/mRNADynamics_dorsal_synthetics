function fluoOverTimePerBin(enhancerName,NotEnhancerName,Color,fiducialTime,numBins)

%% load everything
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

% we'll use the struct called 'combinedCompiledProjects_allEnhancers', wich
% contains one entry per nucleus.
% First, make a smaller struct called 'enhancerStruct'containing only the nc12 nuclei for the enhancer
% specified by the 'DataTye' argument.
for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

%% define rules for what data we'll consider in the analysis and plots
dataTypeIncludeRule = contains({combinedCompiledProjects_allEnhancers.dataSet},enhancerName);
dataTypeExcludeRule = ~contains(lower({combinedCompiledProjects_allEnhancers.dataSet}),lower(NotEnhancerName));
dataTypeRule = dataTypeIncludeRule & dataTypeExcludeRule;

% tempRule = ({combinedCompiledProjects_allEnhancers.dataSet} == "1Dg11_2xDl" |...
%       {combinedCompiledProjects_allEnhancers.dataSet} == "1Dg11_FFF"  );
% dataTypeRule = dataTypeRule & tempRule;

enhancerStruct = combinedCompiledProjects_allEnhancers(dataTypeRule &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

%% bin nuclei 

% ****IMPORTANT!!!**** this line uses the standard dorsal fluorescence, not the
% arbitrary one at a given 'fiducial time'
if isempty(fiducialTime)
    nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];
else
    %this function calculates the Dorsal fluorescence at some arbitrary time
    %in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
    enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);
    nucleiFluorescence = [enhancerStruct.DorsalFluoArbitraryTime];
end

binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end

coveredBins = unique([enhancerStruct.dorsalFluoBin2]);
binValues = binValues(coveredBins); %for plots


%% loop over bins and grab all nuclei
frameRate = 10; % seconds
totalTime = 13; %minutes into NC12
timeArray = 0:10:totalTime*60; % in seconds
Palette = plasma(length(coveredBins));
Palette = cbrewer('seq', 'YlGn', length(coveredBins));


figure
hold on
for b = 1:length(coveredBins)
    numParticles = 0;
    binID = coveredBins(b);
    binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin2]== binID);
    ValidNuclei = length(binStruct);
    AllParticlesFluos = nan(ValidNuclei,length(timeArray));
    % loop over particles to populate ALLPARTICLESFLUOS
    for p = 1:length(binStruct)
        if ~isempty([binStruct(p).particleTimeOn]) %&& ~isempty([binStruct(p).dorsalFluoFeature])
            numParticles = numParticles+1;
            particleTime = binStruct(p).particleTimeSinceAnaphase; %in min
            particleTime = particleTime*60; %in seconds
            particleFluo = binStruct(p).particleFluo;
            particleNCFluo = populateTimeArray(timeArray,particleTime,particleFluo,frameRate);
            AllParticlesFluos(p,:) = particleNCFluo;
        elseif isempty([binStruct(p).particleTimeOn]) %&& ~isempty([binStruct(p).dorsalFluoFeature])
            AllParticlesFluos(p,:) = 0;
        end
    end
    forLegend{b} = num2str(binValues(b));
    shadedErrorBar(timeArray/60,mean(AllParticlesFluos,1),std(AllParticlesFluos,[],1)/sqrt(numParticles),...
        'lineProps',{'Color',Palette(b,:),'LineWidth',2})    
end
hold off
xlabel('time into NC12 (min)')
ylabel('mean spot fluorescence (all nuclei)')
title(enhancerName)
legend(forLegend)
ylim([0 120])


end



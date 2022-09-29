function DorsalFluoTraces = getDorsalFluoTraces(numBins,fiducialTime)

%% load everything
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

% we'll use the struct called 'combinedCompiledProjects_allEnhancers', wich
% contains one entry per nucleus.
% First, make a smaller struct called 'enhancerStruct'containing only the
% nc12 nuclei.

for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

allEnhancersNC12 = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12 &...
    ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));


%% bin nuclei 

% ****IMPORTANT!!!**** this line uses the standard dorsal fluorescence, not the
% arbitrary one at a given 'fiducial time'
if isempty(fiducialTime)
    nucleiFluorescence = [allEnhancersNC12.dorsalFluoFeature];
else
    %this function calculates the Dorsal fluorescence at some arbitrary time
    %in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
    enhancerStruct = DorsalFluoArbitraryTime(allEnhancersNC12,fiducialTime);
    nucleiFluorescence = [allEnhancersNC12.DorsalFluoArbitraryTime];
end

binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(allEnhancersNC12)
    allEnhancersNC12(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end

coveredBins = unique([allEnhancersNC12.dorsalFluoBin2]);
binValues = binValues(coveredBins);

%% bin nuclei in another way
binLimits = linspace(0,4000,numBins);
for n = 1:length(allEnhancersNC12)
    FluoFeature = allEnhancersNC12(n).dorsalFluoFeature;
    bin = find(FluoFeature-binLimits<0,1,'first')-1;
    allEnhancersNC12(n).dorsalFluoBin3 = bin;
end

        
    
%% populate the output struct
ncDuration = 12; %min
frameRate = 10; %seconds
interpPoints = 100;
timeSteps = ceil((ncDuration*60)./frameRate);
numNuclei = 100;
absTime = 0:frameRate:timeSteps*frameRate;
% Palette = cbrewer('seq', 'YlGn', numBins);
Palette = viridis(numBins);

figure
hold on
for bin = 1:numBins-1
    
    n = 1;
    for k = 1:length(allEnhancersNC12)
        if allEnhancersNC12(k).dorsalFluoBin3 == bin
            binStruct(n) =  allEnhancersNC12(k);
            n = n + 1;
        end
    end
    
    numNuclei = min([numNuclei,length(binStruct)]);
    DorsalFluoArray = nan(length(absTime),numNuclei);
    randomlyChosenNuclei = randi([1 length(binStruct)],1,numNuclei);
    
    for n = 1:numNuclei
        
        nucleusTime = binStruct(n).nuclearTimeSinceAnaphase.*60; %it was in min
        nucleusFluo = binStruct(n).dorsalFluoTimeTrace;
%        interpTime = linspace(min(nucleusTime),max(nucleusTime),interpPoints);
%        interpFluo = interp1(nucleusTime,nucleusFluo,interpTime);
        
        for t = 1:length(absTime)
            absoluteTime = absTime(t);
            [delta,nearestTimeIdx] = min(abs(nucleusTime-absoluteTime));
            fluoAtThatTime = nucleusFluo(nearestTimeIdx);
            if delta < 20 %seconds
                DorsalFluoArray(t,n) = fluoAtThatTime;
            end
            
        end
        %[n/numNuclei]
    end 
    
    meanNucleusFluo = nanmean(DorsalFluoArray,2);
    errorNucleusFluo = nanstd(DorsalFluoArray,[],2)./sqrt(numNuclei);
%     errorbar(absTime/60,meanNucleusFluo,errorNucleusFluo,'CapSize',0,'Color',Palette(bin,:),'LineWidth',2)
    
    DorsalFluoTraces(bin).absoluteTime = absTime;
    DorsalFluoTraces(bin).meanDorsalFluo = smooth(meanNucleusFluo);
    DorsalFluoTraces(bin).originalFluoFeature = mean([binStruct.dorsalFluoFeature]);
    DorsalFluoTraces(bin).bin = bin;
    DorsalFluoTraces(bin).binValues = binValues;
    
    clear binStruct
    
end
hold off


        
        
    
    

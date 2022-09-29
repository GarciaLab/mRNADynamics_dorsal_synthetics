function dorsalTimeTraces = generateDorsalTimeTraces(dorsalVals)

% for consistency, use this
% dorsalVals = linspace(0,3800,18)
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
binLimits = dorsalVals;
for n = 1:length(allEnhancersNC12)
    FluoFeature = allEnhancersNC12(n).dorsalFluoFeature;
    bin = find(FluoFeature-binLimits<0,1,'first')-1;
    allEnhancersNC12(n).dorsalFluoBin3 = bin;
end

%% populate the output struct

ncDuration = 12; %min
frameRate = 1; %seconds
interpPoints = 500;
timeSteps = ceil((ncDuration*60)./frameRate);
absTime = 0:frameRate:timeSteps*frameRate;
% Palette = cbrewer('seq', 'YlGn', numBins);
Palette = viridis(length(dorsalVals));

figure
hold on
for bin = 1:(length(dorsalVals)-1)
    
    % check if the nucleus is in this bin and add it to a new struct called
    % binStruct
    n = 1;
    for k = 1:length(allEnhancersNC12)
        if allEnhancersNC12(k).dorsalFluoBin3 == bin
            binStruct(n) =  allEnhancersNC12(k);
            n = n + 1;
        end
    end
    
    % average the time traces of nuclei in this bin
    numNuclei = length(binStruct);
    DorsalFluoArray = nan(length(absTime),numNuclei); % to populate with fluo values and then average over nuclei
    
    for n = 1:numNuclei
        
        nucleusTime = binStruct(n).nuclearTimeSinceAnaphase.*60; %it is originally in min
        nucleusFluo = binStruct(n).dorsalFluoTimeTrace;
        interpTime = linspace(min(nucleusTime),max(nucleusTime),interpPoints);
        interpFluo = interp1(nucleusTime,nucleusFluo,interpTime);

        for t = 1:length(absTime)
            absoluteTime = absTime(t);
            [delta,nearestTimeIdx] = min(abs(interpTime-absoluteTime));
            fluoAtThatTime = interpFluo(nearestTimeIdx);
            DorsalFluoArray(t,n) = fluoAtThatTime;       
        end
        
    end 
    
    meanNucleusFluo = nanmean(DorsalFluoArray,2);
    %errorNucleusFluo = nanstd(DorsalFluoArray,[],2)./sqrt(numNuclei);
    %errorbar(absTime/60,meanNucleusFluo,errorNucleusFluo,'CapSize',0,'Color',Palette(bin,:),'LineWidth',2)
    
    DorsalFluoTraces(bin).absoluteTime = absTime;
    DorsalFluoTraces(bin).meanDorsalFluo = smooth(meanNucleusFluo);
    DorsalFluoTraces(bin).originalFluoFeature = mean([binStruct.dorsalFluoFeature]);
    DorsalFluoTraces(bin).bin = bin;
    DorsalFluoTraces(bin).binValue = (dorsalVals(bin)+dorsalVals(bin+1))/2;
    
    clear binStruct
    
end

Path = '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';
save([Path '/DorsalFluoTraces.mat'],'DorsalFluoTraces')


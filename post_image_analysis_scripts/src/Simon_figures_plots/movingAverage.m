function [embryoMeanOverAll embryoSEMOverAll] = movingAverage(DataType,metric,averagingWindow,fiducialTime,Color,ax)
% DataType should be an enhancer such as '1Dg11'
% metrics can be  maxfluo, accumulatedfluo or fraction
% averaging window is the number of nuclei
% ax is for plotting from the plotEVERYTHING funcion

% set up stuff
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat']);%, 'dorsalResultsDatabase')
% combinedCompiledProjects_allEnhancers is a struct with one nucleus per
% entry for all enhancers.
% AllNC12Struct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);
% AllNC12ThisEnhancer = AllNC12Struct(contains({AllNC12Struct.dataSet}, DataType));
% AllNC12ThisEnhancer = DorsalFluoArbitraryTime(AllNC12ThisEnhancer,fiducialTime);
% 
for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

enhancerStruct = combinedCompiledProjects_allEnhancers(contains({combinedCompiledProjects_allEnhancers.dataSet},DataType)  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));
%this function calculates the Dorsal fluorescence at some arbitrary time
%in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);


% sort the struct according to dorsal fluorescence
tempTable = struct2table(enhancerStruct); % convert the struct to a table
sortedT = sortrows(tempTable, 'dorsalFluoFeature'); % sort the table by dorsal fluo
sortedStruct = table2struct(sortedT); % change it back to struct array 

%get rid of nans
sortedStruct = sortedStruct(~isnan([sortedStruct.dorsalFluoFeature]));


% to store everything
maxFluo = [];
accFluo=[];
fraction=[];
timeon=[];
duration=[];
dorsalFluos = [sortedStruct.dorsalFluoFeature];

% to count the number of spots
NParticles = 0;
    
%loop over nuclei of this enhancer to populate the vectors
for i = 1:length(sortedStruct)    
    isOn(i) = ~isempty(sortedStruct(i).particleAccumulatedFluo);    
    if isOn(i)
        accFluo(i) = sortedStruct(i).particleAccumulatedFluo;
        maxFluo(i) = max(sortedStruct(i).particleFluo95);
        timeOn(i) = sortedStruct(i).particleTimeOn;
        duration(i) = sortedStruct(i).particleDuration;
        NParticles = NParticles+1;
    else
        accFluo(i) = 0;
        maxFluo(i) = 0;
        timeOn(i) = nan;
        duration(i) = 0;
    end

end

% convert 0s to NaNs in maxFluo and accFluo and duration so that we only take active
% nucleiS into account 
accFluo(accFluo==0)=nan;
maxFluo(maxFluo==0)=nan;
duration(duration==0)=nan;

% do a moving window with the same number of nuclei.
dorsalFluoWindow = [];
fractionOnWindow =[];
accFluoWindow = [];
maxFluoWindow = [];
timeOnWindow = [];
durationWindow = [];

for i = 1:length(sortedStruct) 
    windowEnds = min(i+averagingWindow,length(sortedStruct));
    dorsalFluoWindow(i) = mean(dorsalFluos(i:windowEnds));
    fractionOnWindow(i) = mean(isOn(i:windowEnds));
    accFluoWindow(i) = nanmean(accFluo(i:windowEnds));
    maxFluoWindow(i) = nanmean(maxFluo(i:windowEnds));
    timeOnWindow(i) = nanmean(timeOn((i:windowEnds)));
    durationWindow(i) = nanmean(duration((i:windowEnds)));
end

% do a moving window with the same 'jump' in dorsal fluorescence every time
dorsalFluoWindow2 = [];
fractionOnWindow2 =[];
accFluoWindow2 = [];
maxFluoWindow2 = [];
timeOnWindow2 = [];
dorsalFluoBinWidth = 150;

for i = 1:max(dorsalFluos-(averagingWindow/3))
    windowStarts = i;
    windowEnds = i+dorsalFluoBinWidth;
    windowNuclei = logical((dorsalFluos>windowStarts) .* (dorsalFluos<windowEnds));
    dorsalFluoWindow2(i) = mean(dorsalFluos(windowNuclei));
    fractionOnWindow2(i) = mean(isOn(windowNuclei));
    accFluoWindow2(i) = nanmean(accFluo(windowNuclei));
    maxFluoWindow2(i) = nanmean(maxFluo(windowNuclei));
    timeOnWindow2(i) = nanmean(timeOn(windowNuclei));
end


fluoOnNuclei = dorsalFluos(isOn);
fluoOffNuclei = dorsalFluos(~isOn);

%% plot results

if strcmpi(metric,'maxfluo') 
hold on
plot(ax,dorsalFluos,maxFluo,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
plot(ax,dorsalFluoWindow,maxFluoWindow,'Color',Color,'LineWidth',2)
ylim([0 600])
xlim([0 4000])
embryoMeanOverAll = nanmean(maxFluoWindow);
embryoSEMOverAll = nanstd(maxFluoWindow)/sqrt(NParticles);
    
    
elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
hold on
plot(ax,dorsalFluos,accFluo,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
plot(ax,dorsalFluoWindow,accFluoWindow,'Color',Color,'LineWidth',2)
ylim([0 800])    
xlim([0 4000])
embryoMeanOverAll = nanmean(accFluoWindow);
embryoSEMOverAll = nanstd(accFluoWindow)/sqrt(NParticles);
    


elseif contains(lower(metric),'fraction')
yyaxis left
hold on
histogram(fluoOffNuclei,'BinWidth',200,'FaceColor',[.7 .7 .7],'EdgeColor','none')
histogram(fluoOnNuclei,'BinWidth',200,'FaceColor',Color,'EdgeColor','none')
%plot(ax,dorsalFluos,isOn,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
yyaxis right
plot(ax,dorsalFluoWindow,fractionOnWindow,'Color',Color,'LineWidth',2)
ylim([0 1.1])
xlim([0 4000])
embryoMeanOverAll = nanmean(fractionOnWindow);
embryoSEMOverAll = nanstd(fractionOnWindow)/sqrt(sum(isOn));
%plot(dorsalFluoWindow2,fractionOnWindow2,'r')
%hold off

elseif contains(lower(metric),'timeon')
hold on
plot(ax,dorsalFluos,timeOn,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
plot(ax,dorsalFluoWindow,timeOnWindow,'Color',Color,'LineWidth',2) 
ylim([0 10])
xlim([0 4000])
embryoMeanOverAll = nanmean(timeOnWindow);
embryoSEMOverAll = nanstd(timeOnWindow)/sqrt(sum(~isnan(timeOnWindow)));

elseif strcmpi(metric,'duration') 
hold on
plot(ax,dorsalFluos,duration,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
plot(ax,dorsalFluoWindow,durationWindow,'Color',Color,'LineWidth',2)
ylim([0 10])
xlim([0 4000])
embryoMeanOverAll = nanmean(durationWindow);
embryoSEMOverAll = nanstd(durationWindow)/sqrt(sum(~isnan(durationWindow)));
    

end









end



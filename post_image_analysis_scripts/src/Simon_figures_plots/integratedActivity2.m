function integratedActivity2(DataType,metric,fiducialTime,numBins,Color,ax)
% DataType should be an enhancer such as '1Dg11'
% metrics can be  maxfluo, accumulatedfluo or fraction
% ax is for plotting from the plotEVERYTHING funcion

% set up stuff
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat']);%, 'dorsalResultsDatabase')
% combinedCompiledProjects_allEnhancers is a struct with one nucleus per
% entry for all enhancers.
% AllNC12Struct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);
% enhancerStruct = AllNC12Struct(contains({AllNC12Struct.dataSet}, DataType));
for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

enhancerStruct = combinedCompiledProjects_allEnhancers(contains({combinedCompiledProjects_allEnhancers.dataSet},DataType)  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

% bin nuclei 
%this function calculates the Dorsal fluorescence at some arbitrary time
%in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);
nucleiFluorescence = [enhancerStruct.DorsalFluoArbitraryTime];
nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];

binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end


% bin nuclei according to Dorsal-Venus fluorescence
DorsalFluos = [enhancerStruct.dorsalFluoFeature];
binLimits = (0:100:nanmax(DorsalFluos));
BinnedDorsalFluos = BinData(DorsalFluos,binLimits);
BinnedDorsalFluos(isnan(DorsalFluos)) = nan;

for i = 1:length(enhancerStruct)
    enhancerStruct(i).dorsalFluoBin2 = BinnedDorsalFluos(i);
end

% now go over bins and calculate the fraction of active nuclei

FractionPerBin = [];
MaxFluoPerBin = [];
AccumulatedFluoPerBin = [];
TimeOnPerBin = [];
ParticleDurationPerBin = [];

for bin = 1:max(BinnedDorsalFluos)
    BinOnNuclei = 0;
    BinOffNuclei = 0;
    BinMaxFluo = nan;
    BinAccumulatedFluo = nan;
    BinTimeOn = nan;
    BinSpotDuration = nan;
    
    for n = 1:length(enhancerStruct)
        if enhancerStruct(n).dorsalFluoBin2 == bin
            if ~isempty(enhancerStruct(n).particleFluo)
                BinOnNuclei = BinOnNuclei+1;
                BinMaxFluo = [BinMaxFluo enhancerStruct(n).particleFluo95];
                BinAccumulatedFluo = [BinAccumulatedFluo enhancerStruct(n).particleAccumulatedFluo];
                BinTimeOn = [BinTimeOn enhancerStruct(n).particleTimeOn];
                BinSpotDuration = [BinSpotDuration enhancerStruct(n).particleDuration];
            else
                BinOffNuclei = BinOffNuclei+1;
            end
        end
    end
    
    FractionPerBin = [FractionPerBin BinOnNuclei/(BinOnNuclei+BinOffNuclei)];
    MaxFluoPerBin = [MaxFluoPerBin nanmean(BinMaxFluo)];
    AccumulatedFluoPerBin = [AccumulatedFluoPerBin nanmean(BinAccumulatedFluo)];
    TimeOnPerBin = [TimeOnPerBin nanmean(BinTimeOn)];
    ParticleDurationPerBin = [ParticleDurationPerBin nanmean(BinSpotDuration)];
end


integrated_fraction = cumsum(FractionPerBin,'omitnan');
integrated_maxfluo = cumsum(MaxFluoPerBin,'omitnan');
integrated_accumulatedfluo = cumsum(AccumulatedFluoPerBin,'omitnan');
integrated_timeon = cumsum(TimeOnPerBin,'omitnan');
integrated_duration = cumsum(ParticleDurationPerBin,'omitnan');



if contains(lower(metric),'fraction')
plot(binLimits,integrated_fraction,'Color',Color,'LineWidth',2)
ylabel('fraction active')

elseif strcmpi(metric,'maxfluo') 
plot(binLimits,integrated_maxfluo,'Color',Color,'LineWidth',2)
ylabel('maximum fluo')

elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
plot(binLimits,integrated_accumulatedfluo,'Color',Color,'LineWidth',2)
ylabel('accumulated fluo')

elseif contains(lower(metric),'timeon')
plot(binLimits,integrated_timeon,'Color',Color,'LineWidth',2)
ylabel('time on')

elseif contains(lower(metric),'duration')
plot(binLimits,integrated_duration,'Color',Color,'LineWidth',2)
ylabel('duration')


% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % sort the struct according to dorsal fluorescence
% tempTable = struct2table(enhancerStruct); % convert the struct to a table
% sortedT = sortrows(tempTable, 'dorsalFluoFeature'); % sort the table by dorsal fluo
% sortedStruct = table2struct(sortedT); % change it back to struct array 
% 
% %get rid of nans
% sortedStruct = sortedStruct(~isnan([sortedStruct.dorsalFluoFeature]));
% 
% 
% % to store everything
% maxFluo = [];
% accFluo=[];
% fraction=[];
% timeOn = [];
% isOn=[];
% dorsalFluosForOnNuclei=[];
% dorsalFluos = [sortedStruct.dorsalFluoFeature];
% 
%     
% %loop over nuclei of this enhancer to populate the vectors
% counter=1;
% for i = 1:length(sortedStruct)    
%     isOn(i) = ~isempty(sortedStruct(i).particleAccumulatedFluo);    
%     if isOn(i)
%         accFluo(counter) = sortedStruct(i).particleAccumulatedFluo;
%         maxFluo(counter) = max(sortedStruct(i).particleFluo95);
%         timeOn(counter) = sortedStruct(i).particleTimeOn;
%         dorsalFluosForOnNuclei(counter) = dorsalFluos(i);
%         counter = counter+1;
%     end
% end
% 
% % convert 0s to NaNs in maxFluo and accFluo so that we only take active
% % nucleiS into account 
% accFluo(accFluo==0)=nan;
% maxFluo(maxFluo==0)=nan;
% 
% int_maxFluo = cumsum(maxFluo,'omitnan');
% int_accFluo = cumsum(accFluo,'omitnan');
% int_fraction = cumsum(isOn);
% int_timeon = cumsum(timeOn);
% 
% 
% % we want the integration window to be always the same number of nuclei
% 
% %deal with the active nuclei only first
% int_maxFluo2=maxFluo(1);
% int_accFluo2=accFluo(1);
% int_timeon2=timeOn(1);
% dorsalFluosForOnNuclei2=dorsalFluosForOnNuclei(1);
% for j = 0:window:length(maxFluo)-window
%     int_maxFluo2 = [int_maxFluo2 int_maxFluo2(end)+sum(maxFluo(j+1:j+window))];
%     int_accFluo2 = [int_accFluo2 int_accFluo2(end)+sum(accFluo(j+1:j+window))];
%     int_timeon2 = [int_timeon2 int_timeon2(end)+sum(timeOn(j+1:j+window))];
%     dorsalFluosForOnNuclei2 = [dorsalFluosForOnNuclei2 mean(dorsalFluosForOnNuclei(j+1:j+window))];
% end
% %deal with all nuclei now
% int_fraction2=isOn(1);
% dorsalFluosForAll2=dorsalFluos(1);
% for j = 0:window:length(dorsalFluos)-window
%     int_fraction2 = [int_fraction2 int_fraction2(end)+sum(isOn(j+1:j+window))];
%     dorsalFluosForAll2 = [dorsalFluosForAll2 mean(dorsalFluos(j+1:j+window))];
% end
% 
% 
% 
% % plot results
% if strcmpi(metric,'maxfluo') 
% plot(ax,dorsalFluosForOnNuclei2,int_maxFluo2,'b','LineWidth',2)
%     
% elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
% plot(ax,dorsalFluosForOnNuclei2,int_accFluo2,'b''LineWidth',2)
%     
% elseif contains(lower(metric),'fraction')
% plot(ax,dorsalFluosForAll2,int_fraction2,'b','LineWidth',2)
% 
% elseif contains(lower(metric),'timeon')
% plot(ax,dorsalFluosForOnNuclei2,int_timeon2,'b','LineWidth',2)
%     
%     
% end
% 
% 







end





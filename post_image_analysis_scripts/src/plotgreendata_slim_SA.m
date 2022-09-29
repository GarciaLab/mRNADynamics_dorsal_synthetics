function plotgreendata_slim(varargin)

if nargout > 0
    displayFigures = false;
else
    displayFigures = true;
end
    
%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

numBins = 17;

fiducialTime = 6; %mins

for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

%just gonna grab 1dg11_2xDl data
enhancerStruct = combinedCompiledProjects_allEnhancers(...
    [combinedCompiledProjects_allEnhancers.cycle]==12 &...
    ({combinedCompiledProjects_allEnhancers.dataSet} == "1DG_VW_2xDl_FFF" |...
      {combinedCompiledProjects_allEnhancers.dataSet} == "1DG_VW_2xDl_FFF"  )&...
    ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));


% bin nuclei
%this function calculates the Dorsal fluorescence at some arbitrary time
%in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);
nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];

% ****IMPORTANT!!!**** this line below is to use the standard dorsal fluorescence, not the
% arbitrary one at a given 'fiducial time' obtained from DorsalFluoArbitraryTime
nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];

binValues = linspace(0,3789.5,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end
coveredBins = unique([enhancerStruct.dorsalFluoBin2]);
binValues = binValues(coveredBins);

% now make one struct per bin
% define some filters
minEmbryosPerBin = 3;
minNucleiPerEmbryoPerBin = 1;
minOnset = 2; % (min) earliest possible spot detection time to be counted
maxOnset = 8; %(min) latest possible spot detection time to be counted

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

OnNucleiPerBin = [];
OffNucleiPerBin = [];

FractionsPerEmbryo = nan(50,length(coveredBins));
TimeOnsPerEmbryo = nan(50,length(coveredBins));
MaxFluoPerEmbryo = nan(50,length(coveredBins));
AccFluoPerEmbryo = nan(50,length(coveredBins));
DurationPerEmbryo = nan(50,length(coveredBins));


for b = 1:length(coveredBins)
    binID = coveredBins(b);
    binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin2]== binID);
    activeNuc_Bin = length([binStruct.particleTimeOn]);
    inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;
    % take the means and errors across nuclei first because it's easy      
    mean_maxFluo_acrossNuclei_perBin(b) = nanmean([binStruct.particleFluo95]);
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
            maxFluo_perEmbryo(e) = nanmean([embryoStruct.particleFluo95]);
            accFluo_perEmbryo(e) =  nanmean([embryoStruct.particleAccumulatedFluo]);
            fraction_perEmbryo(e) = length([embryoStruct.particleTimeOn])/length(embryoStruct);
            % filter the onset times to 2<true<8 (min since anaphase)
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
    se__totalRNA_acrossEmbryos_perBin(b) = nanstd(totalRNA_perEmbryo)./sqrt(numEmbryos);
    
    %save the single embryo datum for later use
    FractionsPerEmbryo(1:length(fraction_perEmbryo),b) = fraction_perEmbryo;
    TimeOnsPerEmbryo(1:length(timeOn_perEmbryo),b) = timeOn_perEmbryo;
    MaxFluoPerEmbryo(1:length(maxFluo_perEmbryo),b) = maxFluo_perEmbryo;
    AccFluoPerEmbryo(1:length(accFluo_perEmbryo),b) = accFluo_perEmbryo;
    DurationPerEmbryo(1:length(fraction_perEmbryo),b) = duration_perEmbryo;
    %clear  nuclei_per_Embryo maxFluo_perEmbryo accFluo_perEmbryo fraction_perEmbryo timeOn_perEmbryo
     
end

%% scatter of mean onset times (time on) on X vs mean fraction active on Y, per embryo.

dat_fraction_dl =  FractionsPerEmbryo(:);
dat_onset_dl = TimeOnsPerEmbryo(:);
dat_onset_dl_noNan = dat_onset_dl(~isnan(dat_onset_dl(:)));
dat_fraction_dl = dat_fraction_dl(~isnan(dat_onset_dl(:)));
dat_onset_dl = dat_onset_dl_noNan;

c = binValues; 
Palette = brewermap(length(c),'Greens');

figure
hold on
for bin = 1:b
    binFractionData = FractionsPerEmbryo(:,bin);
    binOnsetData = TimeOnsPerEmbryo(:,bin);
    Xdata = binOnsetData(~isnan(binOnsetData));
    Ydata = binFractionData(~isnan(binOnsetData));
    Color = Palette(bin,:);
    plot(binOnsetData,binFractionData,'o','MarkerFaceColor',Color,'MarkerSize',11,'MarkerEdgeColor','k')
end
hold off
xlim([0 10])
ylim([0 1.1])
xlabel('onset time (min)')
ylabel('fraction of active nuclei')

%% scatter of max fluo on X vs mean fraction active on Y, per embryo.
% 
% dat_fraction_dl =  FractionsPerEmbryo(:);
% dat_maxfluo_dl = MaxFluoPerEmbryo(:);
% dat_maxfluo_dl_noNan = dat_maxfluo_dl(~isnan(dat_maxfluo_dl(:)));
% dat_fraction_dl = dat_fraction_dl(~isnan(dat_maxfluo_dl(:)));
% dat_maxfluo_dl = dat_maxfluo_dl_noNan;

c = binValues; 
Palette = brewermap(length(c),'Greens');

figure
hold on
for bin = 1:b
    binFractionData = FractionsPerEmbryo(:,bin);
    binmaxfluoData = MaxFluoPerEmbryo(:,bin);
    Xdata = binmaxfluoData(~isnan(binmaxfluoData));
    Ydata = binFractionData(~isnan(binOnsetData));
    Color = Palette(bin,:);
    plot(binmaxfluoData,binFractionData,'o','MarkerFaceColor',Color,'MarkerSize',11,'MarkerEdgeColor','k')
end
hold off
xlim([0 600])
ylim([0 1.1])
xlabel('max fluo (a.u)')
ylabel('fraction of active nuclei')



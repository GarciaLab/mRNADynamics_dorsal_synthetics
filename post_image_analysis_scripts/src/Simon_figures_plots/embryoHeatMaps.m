function embryoHeatMaps(prefix,numBins)


%% load stuff
FrameInfo = getFrameInfo(LiveExperiment(prefix));

CompiledParticles = getCompiledParticles(LiveExperiment(prefix));
CompiledParticles = CompiledParticles.CompiledParticlesStruct;

Schnitzcells = getSchnitzcells(LiveExperiment(prefix));

for s = 1:length(Schnitzcells)
    Schnitzcells(s).originalSchnitzIdx = s;
end

NC12Schnitzcells = Schnitzcells([Schnitzcells.cycle]==12);
for i = 1:length(NC12Schnitzcells)
    if isempty(NC12Schnitzcells(i).FluoFeature)
        NC12Schnitzcells(i).FluoFeature = nan;
    end
end
NC12Schnitzcells2 = NC12Schnitzcells(~isnan([NC12Schnitzcells.FluoFeature]));

%% get the NC12 frames
NCs = [FrameInfo.nc];
AbsTime = [FrameInfo.Time];
NC12Start = find(NCs==12,1,'first');
NC12End = find(NCs==12,1,'last');
numFrames = NC12End-NC12Start;
AbsNC12Time = AbsTime(NC12Start:NC12End);
AbsNC12Time = AbsNC12Time(2:end)-NC12Start;
%% bin nuclei from scratch
nucleiFluorescence = [NC12Schnitzcells2.FluoFeature];
binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(NC12Schnitzcells2)
    NC12Schnitzcells2(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end
coveredBins = unique([NC12Schnitzcells2.dorsalFluoBin2]);
binValues = binValues(coveredBins);

%% go bin by bin and make a heatmap of the nuclei in that bin

for bi = coveredBins
    figure
    nucleiThisBin = find([NC12Schnitzcells2.dorsalFluoBin2]==bi);
    nucleiPerBin = sum([NC12Schnitzcells2.dorsalFluoBin2]==bi);
    binNucleiIdx = [NC12Schnitzcells2.originalSchnitzIdx];
    binNucleiIdx = binNucleiIdx(nucleiThisBin);
    AllNucleiSpotArray = zeros(numFrames,nucleiPerBin);
    counter=1;
    for n = binNucleiIdx
        particleIdx = find([CompiledParticles.schnitz]==n);
        if particleIdx
            particleFluo = CompiledParticles(particleIdx).Fluo;
            particleFrames = CompiledParticles(particleIdx).Frame;
            particleFrames = particleFrames-NC12Start;
            AllNucleiSpotArray(particleFrames,counter) = particleFluo;
        else
            AllNucleiSpotArray(:,counter) = 0;
        end
        counter = counter+1;
    end
    imagesc(AllNucleiSpotArray')
    title(['bin = ' num2str(bi)])
end


%% pool together all bins (all nucleaises)

%figure
%nucleiThisBin = find([NC12Schnitzcells2.dorsalFluoBin2]==bi);
%nucleiPerBin = sum([NC12Schnitzcells2.dorsalFluoBin2]==bi);
%binNucleiIdx = [NC12Schnitzcells2.originalSchnitzIdx];
%binNucleiIdx = binNucleiIdx(nucleiThisBin);
AllNucleiSpotArray = zeros(length(NC12Schnitzcells2),numFrames);
AllNucleiDorsalFluoVec = zeros(length(NC12Schnitzcells2),1);

for n = 1:length(NC12Schnitzcells2)
    SchnitzID = NC12Schnitzcells2(n).originalSchnitzIdx;
    particleIdx = find([CompiledParticles.schnitz]==SchnitzID);
    AllNucleiDorsalFluoVec(n) = NC12Schnitzcells2(n).FluoFeature;
    if particleIdx
        particleFluo = CompiledParticles(particleIdx).Fluo;
        particleFrames = CompiledParticles(particleIdx).Frame;
        particleFrames = particleFrames-NC12Start;
        AllNucleiSpotArray(n,particleFrames) = particleFluo;
    else
        AllNucleiSpotArray(n,:) = 0;
    end
end
AllNucleiSpotArray(AllNucleiSpotArray<0)=0;

%% sort nuclei by their particle fluorescence
%concatenate the Dorsal fluorescence to the particle fluo array
AllNucleiSpotArray = [AllNucleiSpotArray AllNucleiDorsalFluoVec];
% sume the fluo of each spot over time
sumFluo = sum(AllNucleiSpotArray(:,1:end-1),2);
% concatenate the sum spot fluo in the last column
AllNucleiSpotArray = [AllNucleiSpotArray sumFluo];

% sort array based on sum spot fluo
SortedAllNucleiArray = sortrows(AllNucleiSpotArray,size(AllNucleiSpotArray,2));
SortedAllNucleiArray(:,end)=0;
figure
%subplot(1,2,1)
imagesc(flip(SortedAllNucleiArray(:,1:end-3)))
colorbar
colormap plasma
title(['All nc12 nuclei from ' prefix])
xlabel('time into nc12 (min)')
ylabel('nuclei')

figure
%subplot(1,2,2)
imagesc(flip(SortedAllNucleiArray(:,end-1)))
colorbar
palette = cbrewer('seq', 'YlGn', size(SortedAllNucleiArray,1));
colormap(palette)
title(['Dorsal fluo - all nc12 nuclei from ' prefix])
ylabel('nuclei')
%% deal with labels
% timeLabels = 0:5:15;
% find(AbsNC12Time./60==
% labels = AbsNC12Time(xticks)./60;
% xticklabels(num2str(labels))
%% 








end

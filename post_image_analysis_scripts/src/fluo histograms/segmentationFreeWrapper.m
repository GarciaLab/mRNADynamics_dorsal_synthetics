nullvals = [getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_1'), getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_2'),...
    getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_3'), getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_4')];


wtsets = getProjectPrefixes('1Dg11_2xDl', 'onlyApproved');
wtvals = [];
for k = 1:length(wtsets)
    wtvals = [wtvals, getDoGVals(wtsets{k})]; %#ok<AGROW>
end

segmentationFreeMS2Analysis('', 20, 'vals', nullvals - min(nullvals))
segmentationFreeMS2Analysis(Prefix2, 20,'ax', gca, 'vals', wtvals - min(wtvals))



%%

wtsets = getProjectPrefixes('1dg2x_hist');
wtvals = [];
for k = 1:length(wtsets)
    Spots = getSpots(LiveExperiment(wtsets{k}));
    for frame = 1:length(Spots)
        for s = 1:length(Spots(frame).Fits)
            wtvals = [wtvals, double(Spots(frame).Fits(s).DOGIntensity)]; %#ok<AGROW>
        end
    end
    Spots = [];
end

wtvals = log10((wtvals/100)-100 + 1);


segmentationFreeMS2Analysis('', 20, 'vals', nullvals - min(nullvals))
segmentationFreeMS2Analysis(Prefix2, 20,'ax', gca, 'vals', wtvals - min(wtvals))

%%

wtvals = [];
for k = 1:length(wtsets)
    
    Particles = getParticles(LiveExperiment(wtsets{k}));
    Spots = getSpots(LiveExperiment(wtsets{k}));
    
    for CurrentParticle = 1:length(Particles)
        
        
        for frameIndex=1:length(Particles(CurrentParticle).Frame)
            
            Frame(frameIndex)=Particles(CurrentParticle).Frame(frameIndex);
            spot = Spots(Frame(frameIndex)).Fits(Particles(CurrentParticle).Index(frameIndex));
            zIndex=find(spot.brightestZ == spot.z);
            AmpDogMax(frameIndex) =  double(spot.DOGIntensity(zIndex));
            
        end
        
%         wtvals = [wtvals, AmpDogMax(AmpDogMax > quantile(AmpDogMax(:), .9))];
        wtvals = [wtvals, max(AmpDogMax)];
    end
    
end

wtvals = log10((wtvals/100)-100 + 1);


segmentationFreeMS2Analysis('', 20, 'vals', nullvals - min(nullvals))
segmentationFreeMS2Analysis('', 20,'ax', gca, 'vals', wtvals - min(wtvals))


%%

wtsets = getProjectPrefixes('1Dg11_2xDl', 'onlyApproved');
wtvals = [];
for k = 1:length(wtsets)
    
    liveExperiment = LiveExperiment(wtsets{k});
    Particles = getParticles(liveExperiment);
    Spots = getSpots(liveExperiment);
    
    DogOutputFolder = [liveExperiment.procFolder, 'dogs', filesep];
    dogDir = dir([DogOutputFolder, '*_ch0','*.*']);
    dogDir = {dogDir.name};
    haveStacks = any(cellfun(@(x) ~contains(x, '_z'), dogDir));
    FrameInfo = getFrameInfo(liveExperiment);

    nFrames = length(FrameInfo);

    for CurrentParticle = 1:length(Particles)
        
        
        for frameIndex=1:length(Particles(CurrentParticle).Frame)
            
            x(frameIndex) = Particles(CurrentParticle).xPos(frameIndex);
            y(frameIndex) = Particles(CurrentParticle).yPos(frameIndex);
            z(frameIndex) = Particles(CurrentParticle).zPos(frameIndex);            
            Frame(frameIndex) = Particles(CurrentParticle).Frame(frameIndex);
            
            dog = double(imreadStack([DogOutputFolder, filesep, dogDir{Frame(frameIndex)}]));
            
            ampdog = dog(x(frameIndex), y(frameIndex), z(frameIndex)); 
            wtvals = [wtvals, dog(x(frameIndex), y(frameIndex), z(frameIndex))];
            
        end
        
%         wtvals = [wtvals, AmpDogMax(AmpDogMax > quantile(AmpDogMax(:), .9))];
    end
    
end


wtvals = log10((wtvals/100)-100 + 1);


nullvals = [getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_1', 'nuclearMask'), getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_2','nuclearMask'),...
    getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_3','nuclearMask'), getDoGVals('2019-06-27-1DG_11xDl1MCPmCh_4','nuclearMask')];


out1 = segmentationFreeMS2Analysis('', 20, 'vals', nullvals - min(nullvals));
out2 = segmentationFreeMS2Analysis('', 20,'ax', gca, 'vals', wtvals - min(wtvals));

%find the point of intersection of the two gaussians to
%estimate the fraction of false positives and false negatives
intersection = find(islocalmin(abs(out1.y-out2.y))==1);

plot(gca, out1.xvalues(intersection), out1.y(intersection), 'o');

cdf1 = cdf(out1.pd, out1.xvalues);
cdf2 = cdf(out2.pd, out2.xvalues);

%the false positive rate in percent is aocright*100. this means that % of spots we called spots were
%actually noise. false negative rate is aocleft*100. this is % of spots
%categorized as noise but aren't. 
aocleft = cdf2(intersection);
aocright = cdf1(end) - cdf1(intersection);

FNR = aocleft*100
FPR = aocright*100
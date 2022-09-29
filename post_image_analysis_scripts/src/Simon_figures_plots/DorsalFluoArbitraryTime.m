function struct = DorsalFluoArbitraryTime(struct,fiducialTime)

for n = 1:length(struct)
    timeSinceAnaphase = struct(n).nuclearTimeSinceAnaphase;
    DorsalFluorescence = struct(n).dorsalFluoTimeTrace;
    
    [delta NearestFrameIdx] = min(abs(timeSinceAnaphase-fiducialTime));
    
    if delta < 1 %in minutes   
        struct(n).DorsalFluoArbitraryTime = DorsalFluorescence(NearestFrameIdx);
    else
        struct(n).DorsalFluoArbitraryTime = nan;
    end
end

end 
    
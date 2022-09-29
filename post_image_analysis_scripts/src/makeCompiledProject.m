function compiledProject = makeCompiledProject(Prefix)

[~, resultsFolder] = getDorsalFolders;

liveExperiment = LiveExperiment(Prefix);

load([resultsFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
load([resultsFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'schnitzcells')
load([resultsFolder,filesep,Prefix,filesep,'CompiledParticles.mat'], 'CompiledParticles')

if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{1};
end



if ~isfield(schnitzcells, 'TimeSinceAnaphase')
    
    ncFrames = [zeros(1,8),liveExperiment.anaphaseFrames'];
    schnitzcells = addRelativeTimeToSchnitzcells(...
        schnitzcells, FrameInfo, ncFrames);
    
    CompiledParticles = addRelativeTimeToCompiledParticles(...
        schnitzcells, CompiledParticles, ncFrames, FrameInfo);
end




compiledProject = [];

%do some nuclear quality controle

threshold_frames = 10; %10 frames ~1.5mins under standard conditions
for k = 1:length(schnitzcells)
    if length(schnitzcells(k).frames) <= threshold_frames
        schnitzcells(k).Approved = false;
    else
        schnitzcells(k).Approved = true;
    end
end

approvedSchnitzes = find([schnitzcells.Approved]);

n = 0;
for s = approvedSchnitzes
    
    n = n + 1;
    %nuclear compilation
    compiledProject(n).nuclearFrames = schnitzcells(s).frames;
    compiledProject(n).cycle = schnitzcells(s).cycle;
    if isfield(schnitzcells, 'compiledParticle')
        compiledProject(n).compiledParticle = schnitzcells(s).compiledParticle;
    else
        compiledProject(n).compiledParticle = [];
        schnitzcells(s).compiledParticle = [];
    end
    
    if ~isempty(schnitzcells(s).dorsalFluoBin)
        compiledProject(n).dorsalFluoBin = schnitzcells(s).dorsalFluoBin;
    else
        compiledProject(n).dorsalFluoBin = nan;
    end
    compiledProject(n).dorsalFluoFeature = schnitzcells(s).FluoFeature;
    compiledProject(n).dorsalFluoTimeTrace= schnitzcells(s).FluoTimeTrace;
    compiledProject(n).nuclearTimeSinceAnaphase = schnitzcells(s).timeSinceAnaphase;
    compiledProject(n).prefix = Prefix;
    compiledProject(n).index = n;
    
    
    %particle compilation
    p = schnitzcells(s).compiledParticle;
    if ~isempty(p) && CompiledParticles(p).Approved
        
        
        compiledProject(n).particleFrames = CompiledParticles(p).Frame;
        compiledProject(n).particleTimeSinceAnaphase = CompiledParticles(p).timeSinceAnaphase;
        
        
        
        compiledProject(n).particleFluo= CompiledParticles(p).Fluo;
        compiledProject(n).particleFluo3Slice = CompiledParticles(p).Fluo3;
        compiledProject(n).particleOffset = CompiledParticles(p).Off;
        compiledProject(n).particleFluoError = CompiledParticles(p).FluoError;

        if isfield(CompiledParticles, 'Fluo3DRaw')
            compiledProject(n).particleFluo3DRaw = CompiledParticles(p).Fluo3DRaw;
            compiledProject(n).particleFluo3D = CompiledParticles(p).Fluo3DGauss;

        else
            compiledProject(n).particleFluo3D = nan(1, length(CompiledParticles(p).Frame));            compiledProject(n).particleFluo3D = nan(1, length(CompiledParticles(p).Frame));
            compiledProject(n).particleFluo3DRaw = nan(1, length(CompiledParticles(p).Frame));
            warning("missing gauss 3D intensities for " + Prefix);
        end
        
        %                 fluo = compiledProject(n).particleFluo3D;
        fluo = compiledProject(n).particleFluo3Slice;
        tau = compiledProject(n).particleTimeSinceAnaphase;
        
        %some preliminary additional statistics about the spot
        compiledProject(n).particleDuration = max(tau) - min(tau);
        compiledProject(n).particleFluo95= fluo(fluo>=prctile(fluo,95));
        compiledProject(n).particleTimeOn = min(tau);
        
        if length(fluo) > 1
            compiledProject(n).particleAccumulatedFluo = trapz(tau, fluo, 2);
        else
            compiledProject(n).particleAccumulatedFluo = fluo;
        end
        
        
    else
        
        compiledProject(n).particleFrames = [];
        compiledProject(n).particleFluo= [];
        compiledProject(n).particleFluo3Slice = [];
        compiledProject(n).particleOffset = [];
        compiledProject(n).particleFluo3D = [];
        compiledProject(n).particleFluo3DRaw = [];
        compiledProject(n).particleDuration = [];
        compiledProject(n).particleFluo95= [];
        compiledProject(n).particleTimeOn = [];
        compiledProject(n).particleAccumulatedFluo = [];
        compiledProject(n).particleFluoError = [];
        
    end
end


assert(~isempty(compiledProject))

save([resultsFolder,filesep,Prefix,filesep,'compiledProject.mat'], 'compiledProject');


end
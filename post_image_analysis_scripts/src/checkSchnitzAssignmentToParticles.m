function checkSchnitzAssignmentToParticles(Prefix)
% goes over compiledParticles and checks that the schnitz assigned to each
% particle is the right one


% [~,~,DropboxFolder,~, PreProcPath,...
%     ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
% resultsFolder = [DropboxFolder, filesep, Prefix];

liveExperiment = LiveExperiment(Prefix);
resultsFolder = liveExperiment.resultsFolder;

load([resultsFolder,'CompiledParticles.mat']);
load([resultsFolder, 'FrameInfo.mat']);
load([resultsFolder, Prefix, '_lin.mat'])

channel = 1;
CompiledParticlesStruct = CompiledParticles{channel};

%%
if ~isempty(CompiledParticlesStruct)
    particleNCs = [CompiledParticlesStruct.nc];
    nc12Idx = [1:length(CompiledParticlesStruct)] .* (particleNCs==12);

    minDistance = 50; %pixels of distance between a spot and a nucleus centroid to be matched

    if ~isempty(CompiledParticles)

    for p = 1:nc12Idx

        particleXPosPerFrame = CompiledParticlesStruct(p).xPos;
        particleYPosPerFrame = CompiledParticlesStruct(p).yPos;
        particleFrames = CompiledParticlesStruct(p).Frame;
        closestNucleusPerFrame = nan(1,length(particleFrames));

        for f = 1:length(particleFrames)

            frame = particleFrames(f);
            particleXPosThisFrame = particleXPosPerFrame(f);
            particleYPosThisFrame = particleYPosPerFrame(f);
            distanceToNucleiThisFrame = nan(1,length(schnitzcells));
    %         nucleiThisFrameIDs = [];

            for n = 1:length(schnitzcells)
                schnitzFrames = schnitzcells(n).frames;
                if ismember(frame,schnitzFrames)
                    schnitzXPos = schnitzcells(n).cenx;
                    schnitzXPos = double(schnitzXPos(schnitzFrames==frame));
                    schnitzYPos = schnitzcells(n).ceny;
                    schnitzYPos = double(schnitzYPos(schnitzFrames==frame));
                    distanceToNucleiThisFrame(n) = sqrt((particleXPosThisFrame-schnitzXPos)^2+(particleYPosThisFrame-schnitzYPos)^2);
                end
            end
            [dist,idx] = min(distanceToNucleiThisFrame);
            if dist < minDistance
                closestNucleusPerFrame(f) = idx;
            end
        end
        closestNucleusEver = mode(closestNucleusPerFrame);
        if ~isnan(closestNucleusEver)
            if CompiledParticlesStruct(p).schnitz ~= closestNucleusEver
                disp(['original assigned schnitz is not the right one!' Prefix])
                CompiledParticlesStruct(p).schnitz = closestNucleusEver;
            end   
        end
    end
    end


    %disp('congrats! schnitzes were correctly assigned to particles')

    clear closestNucleusEver closestNucleusPerFrame idx dist minDistance distanceToNucleiThisFrame schnitzYPos schnitzXPos particleYPosThisFrame particleXPosThisFrame frame f n particleNCs nc12Idx

    save([resultsFolder,'CompiledParticles.mat']);
end

%%

            

end
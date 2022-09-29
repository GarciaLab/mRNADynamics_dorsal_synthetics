function editCompiledProject(DataType)
%this is a more lightweight version of checkparticletracking.
%right now the main purpose is to categorize particles into
%basal vs. trapezoidal.
%the function loops over prefixes in a data status sheet

cleanupObj = onCleanup(@myCleanupFun);


thisProject = LiveProject(DataType); %#ok<NOPRT>    
[~, resultsFolder] = getDorsalFolders;

prefixes = thisProject.includedExperimentNames;

for i = 1:length(prefixes)
    
    Prefix = prefixes{i};
    load([resultsFolder,filesep,Prefix,filesep,'compiledProject.mat'], 'compiledProject');
    
    %we're gonna loop over a subset of this structure, so
    %the index will help us keep track of the identity of our
    %particles in the original struct
    if ~isfield(compiledProject, 'index')
        for k = 1:length(compiledProject)
            compiledProject(k).index = k;
        end
    end
    
    %we just want to look at particles, so let's make that struct here.
    particleIndex = arrayfun(@(x) ~isempty(x.particleFrames),compiledProject);
    particles = compiledProject(particleIndex);
    
    %for now let's just do nc12
    nc = 12;
    particles = particles([particles.cycle] == nc); 
    
    
    mainFigure = figure;
    mainAxes = axes(mainFigure); %#ok<LAXES>
    xlabel(mainAxes, 'time since anaphase (min)')
    ylabel(mainAxes, '3 slice spot intensity')
    
    currentCharacter = 1;
    currentParticle = 1;
    
    
    while (currentCharacter ~= 'x') && ~isempty(particles)
        
        currentParticleIndex = particles(currentParticle).index;
        
        plot(mainAxes, particles(currentParticle).particleTimeSinceAnaphase,...
            particles(currentParticle).particleFluo3Slice);
        
        if isfield(compiledProject, 'trapezoidStatus') &&...
                ~isempty(compiledProject(currentParticleIndex).trapezoidStatus)
            if strcmpi(compiledProject(currentParticleIndex).trapezoidStatus, 'trapezoid')
                set(gcf,'color','g');
            elseif strcmpi(compiledProject(currentParticleIndex).trapezoidStatus, 'basal')
                set(gcf,'color','r');
            elseif strcmpi(compiledProject(currentParticleIndex).trapezoidStatus, 'questionable')
                set(gcf,'color','y');
                %         else
                %             set(gcf,'color','none');
            end
        else
            set(gcf,'color','w');
        end
        
        title({Prefix;...
            ['nc ', num2str(compiledProject(currentParticleIndex).cycle)];...
            ['Particle ', num2str(currentParticleIndex)]});
        
        waitforbuttonpress;
        currentCharacter=get(mainFigure,'currentcharacter');
        
        
        
        if currentCharacter == 'd' && currentParticle < length(particles)
            currentParticle = currentParticle + 1;
        elseif currentCharacter == 'a' && currentParticle > 1
            currentParticle = currentParticle - 1;
        elseif currentCharacter == 'w'
            compiledProject(currentParticleIndex).trapezoidStatus = 'trapezoid';
        elseif currentCharacter == 's'
            compiledProject(currentParticleIndex).trapezoidStatus = 'basal';
        elseif currentCharacter == 'e'
            compiledProject(currentParticleIndex).trapezoidStatus = 'questionable';
        elseif currentCharacter == 'k'
            keyboard;
        end
        
        
    end
    
    answer = inputdlg('Save and exit?');
    unanswered = true;
    while unanswered
        if startsWith(answer, 'y', 'IgnoreCase', true)
            save([resultsFolder,filesep,Prefix,filesep,'compiledProject.mat'], 'compiledProject');
            save([resultsFolder,filesep,Prefix,filesep,'compiledProject.mat'], 'compiledProject');
            unanswered = false;
        elseif startsWith(answer, 'n', 'IgnoreCase', true)
            unanswered = false;
        end
    end
    
end

end
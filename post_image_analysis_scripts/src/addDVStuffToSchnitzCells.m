function [dlfluobins, dlfluobincounts] =...
    addDVStuffToSchnitzCells(DataType, varargin)
%%
displayFigures = false;
saveFigures = false;

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'displayFigures')
        displayFigures = true;
    elseif strcmpi(varargin{k}, 'saveFigures')
        saveFigures = true;
    end
end

[allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType, 'noCompiledNuclei');


if displayFigures
    tileFig = figure();
    tiledlayout('flow');
    holdFig = figure();
end



ch = 1;

for e = 1:length(allData)
        
    liveExperiment = LiveExperiment(Prefixes{e});

    schnitzcells = getSchnitzcells(liveExperiment); 
    
    CompiledParticles = allData(e).Particles.CompiledParticles;    
    
    %to speed this up, remove unnecessary schnitzcells fields
    schnitzcells = removeSchnitzcellsFields(schnitzcells);
    
    % SA: this only checks that the sizes match I think
    checkSchnitzcellsCompiledParticlesConsistency(...
    schnitzcells,...
    CompiledParticles)

    
    %clear out any existing compiledparticle entries
    if isfield(schnitzcells, 'compiledParticle')
        schnitzcells = rmfield(schnitzcells, 'compiledParticle');
    end
    
    % SA: this gets the schnitz ID of each compiled particle and uses it to
    % add the particle info to the corresponding schnitz. Also adds the DV
    % bin info from the schnitz to the particle.
    for p = 1:length(CompiledParticles{ch})
        
        schnitzInd = CompiledParticles{ch}(p).schnitz;
        
        %not sure where the misassignment to compileparticles happens       
        assert(schnitzInd <= length(schnitzcells));        
        schnitzcells(schnitzInd).compiledParticle = uint16(p);
        
        if isfield(CompiledParticles{ch}(p), 'dvbin')
            schnitzcells(schnitzInd).dvbin = ...
                uint8(CompiledParticles{ch}(p).dvbin);
        end
    end
    
    save([resultsFolder,filesep,Prefixes{e},filesep,Prefixes{e},'_lin.mat'], 'schnitzcells');
    
    
    schnitzcells = addDlFluoToSchnitzcells(Prefixes{e}); 
    
    % SA: this only checks that the sizes match I think
    checkSchnitzcellsCompiledParticlesConsistency(...
    schnitzcells,...
    CompiledParticles)

    if displayFigures && saveFigures
        mkdir([resultsFolder, filesep, DataType]);
        
        saveas(tileFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTracesTile.png']);
        saveas(tileFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTracesTile.fig']);
        saveas(tileFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTracesTile.eps']);
        
        saveas(holdFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTraces.png']);
        saveas(holdFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTraces.fig']);
        saveas(holdFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTraces.eps']);
    end
    
   
    
end


end

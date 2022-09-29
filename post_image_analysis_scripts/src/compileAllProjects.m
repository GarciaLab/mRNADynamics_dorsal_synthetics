function dorsalResults = compileAllProjects(DataType)
% 
% dataSets = { '1Dg11_2xDl', '1DgW_2x_Leica',...
%      '1Dg-5_2xDl',...
%     '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl', '1DgSVW2_2xDl','1DgVW_2xDl',...
%     '1DgW_FFF', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF','1Dg-8D_FFF'};
% 
% for k = 1:length(dataSets)
%     compileAllProjects(dataSets{k})
%     
% end

dataTypes = {'1Dg-8D_FFF', '1Dg11_2xDl', '1DgW_2x_Leica',...
    '1DgW_FFF', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF', '1Dg-5_2xDl',...
    '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl', '1DgSVW2_2xDl','1DgVW_2xDl','2Dgc_FFF_2x','TwiPEv5(7)_2xDl'};


thisProject = LiveProject(DataType) %#ok<NOPRT>

[~, resultsFolder] = getDorsalFolders;

% dataTypes = {'1Dg-8D_FFF', '1Dg11_2xDl', '1DgW_2x_Leica',...
%     '1DgW_FFF', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF', '1Dg-5_2xDl',...
%     '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl', '1DgSVW2_2xDl','1DgVW_2xDl'};



prefixes = thisProject.includedExperimentNames;

compiledProjects = cell(1, length(prefixes));
% 
% 
if true
for k = 1:length(prefixes)
%     prefixes{k}
%     clear LiveExperiment
%     TrackNuclei(prefixes{k},'retrack', 'nWorkers', 1);
     integrateSchnitzFluo(prefixes{k});
% %     TrackmRNADynamics(prefixes{k});
%    CompileParticles(prefixes{k},  'minBinSize', 0, 'MinParticles', 0,...
%        'yToManualAlignmentPrompt');
%     fit3DGaussiansToAllSpots(prefixes{k}, 1, 'nWorkers', 1)
%     CompileParticles(prefixes{k},  'minBinSize', 0, 'MinParticles', 0,...
%         'yToManualAlignmentPrompt');
%    alignCompiledParticlesByAnaphase(prefixes{k});
end
end

%% Validate experiments included in the analysis
hasAllPushed = [thisProject.hasSpots, thisProject.hasParticles,...
    thisProject.hasSchnitzcells, thisProject.hasCompiledParticles,...
    thisProject.anaphaseFramesAnnotated];

% assert( all(hasAllPushed) );

% SA I created this function that checks that the schnitz field in
% CompiledParticles contains the right schnitz.
for k = 1:length(prefixes)
    checkSchnitzAssignmentToParticles(prefixes{k})
end
% 
addDVStuffToSchnitzCells(DataType)

binDorsal(DataType, false)


for k = 1:length(prefixes)
    
    compiledProjects{k} = makeCompiledProject(prefixes{k});
    
    if k == 1
        combinedCompiledProjects = compiledProjects{k}; 
    else
        try
            combinedCompiledProjects = [combinedCompiledProjects, compiledProjects{k}]; %#ok<AGROW>
        catch
            [combinedCompiledProjects, compiledProjects{k}] = addFields(combinedCompiledProjects,compiledProjects{k});
            combinedCompiledProjects = [combinedCompiledProjects, compiledProjects{k}]; %#ok<AGROW>
        end
        
    end
    
end

[combinedCompiledProjects.dataSet] = deal(DataType);

save([resultsFolder,filesep,DataType,filesep,'combinedCompiledProjects.mat'], 'combinedCompiledProjects');

% averagedTimeTraces = averageCombinedCompiledProjects(DataType, true);

try
dorsalResults = createDorsalResults(DataType); 
end
% plotDorsalResultsLoop(DataType, 'frac', 1:6, 'hill')
% plotDorsalResultsLoop(DataType, activity)

end
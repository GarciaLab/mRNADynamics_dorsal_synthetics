% dataSets = { '1Dg11_2xDl', '1DgW_2x_Leica',...
%      '1Dg-5_2xDl',...
%     '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl', '1DgSVW2_2xDl','1DgVW_2xDl',...
%     '1DgW_FFF', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF','1Dg-8D_FFF'};
% 

DataType = '1DgW_2x_Leica';
binDorsal(DataType, false, 'dolog')

thisProject = LiveProject(DataType); 
[~, resultsFolder] = getDorsalFolders;
prefixes = thisProject.includedExperimentNames;

for k = 1:length(prefixes) 
   makeCompiledProject(prefixes{k});
end

% createDorsalResults2(DataType, 'minNuclei', 1, 'minEmbryos', 1)
createDorsalResults2(DataType,'minNuclei', 3, 'minEmbryos', 3)

%%
binDorsal(DataType, false)

for k = 1:length(prefixes) 
   makeCompiledProject(prefixes{k});
end

% createDorsalResults2(DataType, 'minNuclei', 1, 'minEmbryos', 1)
createDorsalResults2(DataType,'minNuclei', 3, 'minEmbryos', 3)
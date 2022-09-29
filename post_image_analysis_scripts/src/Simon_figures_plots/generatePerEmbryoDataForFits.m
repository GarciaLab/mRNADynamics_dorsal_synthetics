function [FractionsPerEmbryoAll,OnsetsPerEmbryoAll,enhancerNames,binLimits] = generatePerEmbryoDataForFits(numBins)

%numBins = 18
dorsalVals = linspace(0,3800,numBins);
binWidth = diff(dorsalVals); binWidth = binWidth(1);
binMidValues = dorsalVals(1:end-1) + binWidth/2;

% AffinityDataTypes = unique({'1Dg11_2xDl_FFF','1Dg11_2xDl','1DgW_2x_Leica','1DgW_FFF',...
% '1DgVW_FFF', '1Dg11_FFF','1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1DgSVW2_2xDl',...
% '1DgVW_2xDl','1DgW_2xDl_FFF','1DG_VW_2xDl_FFF'});

enhancerNames = {'1Dg11','1DgS2','1DgW_2x','1DgAW3','1DgSVW2','1DgVVW3','1DgVW'};
excludeString = 'export';


for i = 1:length(enhancerNames)
    dataType = enhancerNames{i};
    [FractionsPerEmbryo,TimeOnsPerEmbryo,MaxFluoPerEmbryo] =  getPerEmbryoDataForFits(dataType,excludeString,dorsalVals);
    FractionsPerEmbryoAll{i} = FractionsPerEmbryo;
    OnsetsPerEmbryoAll{i} = TimeOnsPerEmbryo;
    MaxFluosPerEmbryoAll{i} = MaxFluoPerEmbryo;
end

outputFolder = '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';

save([outputFolder '/FractionsPerEmbryoAll.mat'],'FractionsPerEmbryoAll')   % save variable in the output.mat file
save([outputFolder '/OnsetsPerEmbryoAll.mat'],'OnsetsPerEmbryoAll')   % save variable in the output.mat file
save([outputFolder '/MaxFluoPerEmbryoAll.mat'],'MaxFluosPerEmbryoAll')   % save variable in the output.mat file
save([outputFolder '/enhancerNames.mat'],'enhancerNames')   % save variable in the output.mat file
save([outputFolder '/binMidValues.mat'],'binMidValues')   % save variable in the output.mat file

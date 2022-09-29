
[~, dropboxfolder] = getDorsalFolders;

%1Dg data (highest affinity)
datPath = dropboxfolder + "\manuscript\window\basic\dataForFitting\archive\";
load(datPath + "binMidValues.mat", "binMidValues");
load(datPath + "FractionsPerEmbryoAll.mat", "FractionsPerEmbryoAll");
load(datPath + "OnsetsPerEmbryoAll.mat", "OnsetsPerEmbryoAll");

batchedAffinities = false;


if batchedAffinities
    
    %this setting being true wouldn't make sense, so let's fix it.
    fixKD = false;
    
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
    data_batch = {};
    
    
    for k = 1:numel(enhancers)
        %needs to be nobs x ny
        FractionsPerEmbryo = FractionsPerEmbryoAll{k};
        OnsetsPerEmbryo = OnsetsPerEmbryoAll{k};
        
        
        X = repmat(binMidValues, 1, max(size(FractionsPerEmbryo)));
        F = FractionsPerEmbryo(:)';
        T = OnsetsPerEmbryo(:)';
        X(isnan(F)) = [];
        T(isnan(F)) = [];
        F(isnan(F)) = [];
        data{k}.ydata = [X; F; T]';
    end
    
else
    
    FractionsPerEmbryo = FractionsPerEmbryoAll{1};
    OnsetsPerEmbryo = OnsetsPerEmbryoAll{1};
    %needs to be nobs x ny
    X = repmat(binMidValues, 1, max(size(FractionsPerEmbryo)));
    F = FractionsPerEmbryo(:)';
    T = OnsetsPerEmbryo(:)';
    X(isnan(F)) = [];
    T(isnan(F)) = [];
    F(isnan(F)) = [];
    
    X(isnan(T)) = [];
    F(isnan(T)) = [];
    T(isnan(T)) = [];
    
    data.ydata = [X; F; T]';
end

fraction = F';
onset = T';
x = X';
writetable(table(x, fraction, onset), "C:\Users\owner\Documents\Dorsal-Synthetics-Analysis\dat\"+'basic.csv');

binMidValues = binMidValues';
frac_mean = nanmean(FractionsPerEmbryo, 1)';
frac_onset = nanmean(OnsetsPerEmbryo, 1)';
nFVals = sqrt(length(FractionsPerEmbryo(~isnan(FractionsPerEmbryo))));
nOVals = sqrt(length(OnsetsPerEmbryo(~isnan(OnsetsPerEmbryo))));
frac_ste = (nanstd(OnsetsPerEmbryo, 1)/nOVals)';
onset_ste = (nanstd(FractionsPerEmbryo, 1)/nFVals)';

writetable(table(binMidValues, frac_mean, frac_onset, frac_ste, onset_ste),...
    "C:\Users\owner\Documents\Dorsal-Synthetics-Analysis\dat\"+'binned.csv');




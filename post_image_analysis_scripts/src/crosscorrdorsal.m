close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])
% 
a = combinedCompiledProjects_allEnhancers(strcmpi({combinedCompiledProjects_allEnhancers.dataSet}, '1Dg11_2xDl' )  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12);
% 
% a = combinedCompiledProjects_allEnhancers(...
%     [combinedCompiledProjects_allEnhancers.cycle]==12);

dlfluobins = 0:250:4500;

% b = a(cellfun(@any, {a.particleFrames}));
fluo95 = [];
fluo95_se = [];
for k = 1:length(dlfluobins)
    b = a([a.dorsalFluoBin]==k);
    fluo95(k) = nanmean([b.particleFluo95]);
    fluo95_se(k) = nanstd([b.particleFluo95])./sqrt(length([b.particleFluo95]));
    
end

plot(dlfluobins, fluo95);
errorbar(dlfluobins, fluo95, fluo95_se);

b = a(cellfun(@any, {a.particleFrames}));

for k = 1:length(b)
    nuclearTime = b(k).nuclearTimeSinceAnaphase;
    particleTime = b(k).particleTimeSinceAnaphase;
    f0 = find(particleTime(1)==nuclearTime);
    fend = find(particleTime(end)==nuclearTime);
    dlfluo = [b(k).dorsalFluoTimeTrace];
    dlfluo = dlfluo(f0:fend);
    particleFluo = [b(k).particleFluo3Slice];
    [c,lags]  = xcorr(dlfluo,particleFluo, 'normalized');
    plot(lags, c)
end

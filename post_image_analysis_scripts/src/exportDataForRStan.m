[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
enhancers_1dg = {'1Dg11'};
cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers_1dg{1});
x = dorsalResultsDatabase.dorsalFluoBins(cond);
Y = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
x(isnan(Y)) = [];
Y(isnan(Y)) = [];

writetable(table(x, Y), 'test.csv');

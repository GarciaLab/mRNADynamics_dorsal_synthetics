function plotDorsalPatserScores

patserResultsPath = 'S:\Simon\Dropbox\patser\DV_enhancers';
enhancerFolders = dir(patserResultsPath);

AllScores = nan(100,length(enhancerFolders));
EnhancerNames = {};
count = 1;
for e = 1:length(enhancerFolders)
    folderName = enhancerFolders(e).name;
    if strlength(folderName)>2
        EnhancerNames{count} = folderName;
        resultsPath = [patserResultsPath '\' folderName];
        filePatternMat = fullfile(resultsPath, '*.mat');
        theMatFiles = dir(filePatternMat);
        filePatternTxt = fullfile(resultsPath, '*.txt');
        theTxtFiles = dir(filePatternTxt);
        mapPatserResults([theTxtFiles(1).folder '\' theTxtFiles(1).name]);

        if ~isempty(theMatFiles)
            load([theMatFiles.folder '\' theMatFiles.name],'SiteResults')
            dlScores = [SiteResults.score];
            AllScores(1:length(dlScores),e) = dlScores;
        end
        count=count+1;
    end      
end


figure
hold on
for e = 1:length(enhancerFolders)
    histogram(AllScores(:,e),'DisplayStyle','stairs')
end
hold off


figure
hold on
for e = 1:length(enhancerFolders)
    histogram(AllScores(:,e),'FaceAlpha',0.2,'normalization','probability')
end
hold off

figure
hold on
for e = 1:count
    plot(e,AllScores(:,e),'o')
end
hold off
xticks([1:count])
xticklabels(EnhancerNames)
xtickangle(45)


cutoff = 3.5;
figure
hold on
for e = 1:length(EnhancerNames)
    Scores = AllScores(:,e+2);
    meanScore = nanmean(Scores(Scores>cutoff));
    semScore = nanstd(Scores(Scores>cutoff))./nansum(Scores>cutoff);
    plot(e,Scores(Scores>cutoff),'ro')
    errorbar(e,meanScore,semScore,'ko')
end
hold off
xticks([1:e])
xticklabels(EnhancerNames)
xtickangle(45)
set(gca,'TickLabelInterpreter','none')
ylabel('dl site patser score')


end



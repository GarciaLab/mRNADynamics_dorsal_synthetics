function plotOnsetHistogram()

[~, dropboxfolder] = getDorsalFolders;

%1Dg data (highest affinity)
datPath = dropboxfolder + "\manuscript\window\basic\dataForFitting\archive\";
load(datPath + "OnsetsPerEmbryoAll.mat", "OnsetsPerEmbryoAll");

observations = [];
for k = 1:length(OnsetsPerEmbryoAll)
    observations = [observations; OnsetsPerEmbryoAll{k}(:)];
end
% 
% %use this line for pdf
% figure();
% histogram(observations(:));
% ylabel('counts')
% xlabel('transcription onset time (min)')
% 
% figure();
% histogram(observations(:), 'Normalization', 'cumcount');
% ylabel('cumulative counts')
% xlabel('transcription onset time (min)')
% 

figure();
cdfplot(observations(:));
ylabel('cumulative counts')
xlabel('transcription onset time (min)')
xlim([0, 12])
ylim([0, 1.1])
legend("N="+size(observations(~isnan(observations))))

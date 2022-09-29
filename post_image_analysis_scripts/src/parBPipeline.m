%%
try
    %this has a function that overloads "mad" used in statistics toolbox
    rmpath('.\lib\dependencies\mcmcstat-master\')
    rmpath('P:\Armando\Dorsal-Synthetics-Analysis\lib\mcmcstat-master')
end

%these were segmented using the ParB2-GFP channel
prefix1s_full ={ '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_1','2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_2',...
'2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_3','2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_4',...
'2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_5','2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_6',...
'2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_8','2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_9',...
'2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_10','2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_11',...
'2020-12-05-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_16',...
'2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_18','2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_19',...
'2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_20','2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_21',...
'2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_22','2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_26',...
'2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_27','2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_28',...
'2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_29','2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_30'};

prefix1s ={ '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_5','2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_6',...
'2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_8', '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_10',...
'2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_18','2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_19',...
'2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_20','2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_26',...
'2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_29','2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_30'};

badprefixes = {'2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_1','2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_2',...
'2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_3',...
'2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_27',...
'2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_28',...
'2020-12-05-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_16'...
};


%%
% parpool(numel(prefix1s))
parfor k = 4:length(prefix1s)
    prefix1 = prefix1s{k};
    addSpot2(prefix1);
end
%%
for k = 4:length(prefix1s)
    prefix1 = prefix1s{k};
    fit3DGaussiansToAllSpots(prefix1, 'nWorkers', 20)
end
%%
for k = 1:3
    prefix1 = prefix1s{k};
    segmentSpots(prefix1, [ 6000, 6000], 'segmentChannel', 1, 'nWorkers', 20, 'Shadows', 1, 'nuclearMask', true)
end

%%
% create arrays to store stuff
N = length(prefix1s);
max_intensities_non = cell(N, 1); % for the undetected spots
max_intensities_transcription = cell(N, 1); %detected spots
mean_intensities_non = cell(N, 1); 
mean_intensities_transcription = cell(N, 1);
% loop over datasets and grab the intensities using MeanMax2Channel
parfor k = 4:length(prefix1s)
    
    prefix1 = prefix1s{k};
    
    prefix2 = strcat(prefix1, 'ch2filt');%'ch2filt' is the version of the movie segmented using MCP-mCherry

    
    [max_intensities_non{k},max_intensities_transcription{k},mean_intensities_non{k}, ...
        mean_intensities_transcription{k}] = MeanMax2Channel(prefix1,prefix2);
    
end
%%
 % gather the data and add an offset
MoviesWithMS2 = find(~cellfun(@isempty,max_intensities_transcription));
max_intensities_transcription2 = max_intensities_transcription(MoviesWithMS2);
New_ParB_1Dg_max_trans = [max_intensities_transcription2{:}];
New_ParB_1Dg_max_trans_50 = log(New_ParB_1Dg_max_trans + 200);

max_intensities_non2 = max_intensities_non(MoviesWithMS2);
% max_intensities_non = max_intensities_non(MoviesWithMS2);
% % randomly pick undetected intensities to have the same number as detected ones
% Ndetected = length(New_ParB_1Dg_max_trans_50);
% NpickPerMovie = floor(Ndetected/length(max_intensities_non2));
% for i = 1:length(max_intensities_non2)
%     thisMovieValues = max_intensities_non2{i};
%     pickedValues = randsample(thisMovieValues,NpickPerMovie);
%     max_intensities_non3{i}=pickedValues;
% end 

max_intensities_non3 = max_intensities_non2;

New_ParB_1Dg_max_non = [max_intensities_non3{:}]; 
New_ParB_1Dg_max_non_50 = log(New_ParB_1Dg_max_non + 200);

all_spots = [New_ParB_1Dg_max_non_50, New_ParB_1Dg_max_trans_50];


figure;
histogram(all_spots, 'BinWidth', .06, 'Normalization', 'Probability')
% histogram(all_spots)
set(gca, 'YScale', 'log')

%%
% randomly pick undetected intensities to have the same number as detected
% ones
% Ndetected = length(New_ParB_1Dg_max_trans_50);
% NpickPerMovie = ceil(Ndetected/length(max_intensities_non2));
% for i = 1:length(max_intensities_non2)
%     thisMovieValues = max_intensities_non2{i};
%     pickedValues = randsample(thisMovieValues,NpickPerMovie);
%     max_intensities_non3{i}=pickedValues;
% end

    
    
    





% plot results
figure;
% yyaxis left
h1 = histogram( New_ParB_1Dg_max_trans_50,'EdgeColor','none','FaceColor','b','FaceAlpha',0.7, 'BinWidth', .1, 'Normalization', 'probability');
% h1 = histogram( New_ParB_1Dg_max_trans_50,'FaceColor','none','EdgeColor','b', 'LineWidth', 3, 'DisplayStyle', 'stairs', 'BinWidth', .07);
pd1 = fitdist(New_ParB_1Dg_max_trans_50','normal');
x1 = min(New_ParB_1Dg_max_trans_50):.01:max(New_ParB_1Dg_max_trans_50);
y1 = pdf(pd1,x1');
% plot(x1,y1, 'LineWidth', 3)
% set(gca,'yscale','log')
% xlim([4.07, 6.7])

% h1.BinWidth = 40;

%
hold on
h2 = histogram(New_ParB_1Dg_max_non_50,'FaceColor','#FB9D2D','EdgeColor','none','FaceAlpha',0.7, 'BinWidth', .1, 'Normalization', 'probability');
% h2 = histogram( New_ParB_1Dg_max_non_50,'FaceColor','none','EdgeColor','#FB9D2D', 'LineWidth', 3, 'DisplayStyle', 'stairs', 'BinWidth', .07);

pd2 = fitdist(New_ParB_1Dg_max_non_50','normal');
x2 = min(New_ParB_1Dg_max_non_50):.01:max(New_ParB_1Dg_max_non_50);
y2 = pdf(pd2,x2');
% plot(x2,y2, 'LineWidth', 3)
% set(gca,'yscale','log')
% xlim([4.07, 6.7])
% h2.BinWidth = 40;

% set(gca,'yscale','log')

gm = gmdistribution([pd1.mu, pd2.mu]', permute([.26, .04]', [3 2 1]), [0.393939, 0.606061]')
% gm = gmdistribution([pd1.mu, pd2.mu]', permute([pd1.sigma, pd2.sigma]', [3 2 1]), [0.393939, 0.606061]');


hold on
% yyaxis right
% figure;
% h3 = histogram(all_spots,'BinWidth', .07, 'Normalization', 'probability','EdgeColor','none','FaceColor','b','FaceAlpha',0.7);
% h3 = histogram(all_spots,'BinWidth', .07, 'EdgeColor','none','FaceColor','b','FaceAlpha',0.7);
h3 = histogram(all_spots,'BinWidth', .06, 'FaceColor','none','EdgeColor','k', 'LineWidth', 4, 'DisplayStyle', 'stairs', 'Normalization', 'probability')

% set(gca,'yscale','log')
% figure;
% histfit(all_spots, 40, 'normal')
% pd = fitdist(all_spots','Kernel','Kernel','normal', 'Support', 'positive', 'Width', .015);
pd = fitgmdist(double(all_spots'),2);
hold on
x_values = min(all_spots):.01:max(all_spots);
y = pdf(pd,x_values');
% plot(x_values,y, 'LineWidth', 3)
set(gca,'yscale','log')
% xlim([4.07, 6.7])

%%


% plot results
figure;
% yyaxis left
h1 = histogram( New_ParB_1Dg_max_trans_50,'EdgeColor','none','FaceColor','b','FaceAlpha',0.7, 'BinWidth', .1);
% h1 = histogram( New_ParB_1Dg_max_trans_50,'FaceColor','none','EdgeColor','b', 'LineWidth', 3, 'DisplayStyle', 'stairs', 'BinWidth', .07);
pd1 = fitdist(New_ParB_1Dg_max_trans_50','normal');
x1 = min(New_ParB_1Dg_max_trans_50):.01:max(New_ParB_1Dg_max_trans_50);
y1 = pdf(pd1,x1');
% plot(x1,y1, 'LineWidth', 3)
% set(gca,'yscale','log')
% xlim([4.07, 6.7])

% h1.BinWidth = 40;

%
hold on
h2 = histogram(New_ParB_1Dg_max_non_50,'FaceColor','#FB9D2D','EdgeColor','none','FaceAlpha',0.7, 'BinWidth', .1);
% h2 = histogram( New_ParB_1Dg_max_non_50,'FaceColor','none','EdgeColor','#FB9D2D', 'LineWidth', 3, 'DisplayStyle', 'stairs', 'BinWidth', .07);

pd2 = fitdist(New_ParB_1Dg_max_non_50','normal');
x2 = min(New_ParB_1Dg_max_non_50):.01:max(New_ParB_1Dg_max_non_50);
y2 = pdf(pd2,x2');
% plot(x2,y2, 'LineWidth', 3)
% set(gca,'yscale','log')
% xlim([4.07, 6.7])
% h2.BinWidth = 40;

% set(gca,'yscale','log')

gm = gmdistribution([pd1.mu, pd2.mu]', permute([.26, .04]', [3 2 1]), [0.393939, 0.606061]')
% gm = gmdistribution([pd1.mu, pd2.mu]', permute([pd1.sigma, pd2.sigma]', [3 2 1]), [0.393939, 0.606061]');


hold on
% yyaxis right
% figure;
% h3 = histogram(all_spots,'BinWidth', .07, 'Normalization', 'probability','EdgeColor','none','FaceColor','b','FaceAlpha',0.7);
% h3 = histogram(all_spots,'BinWidth', .07, 'EdgeColor','none','FaceColor','b','FaceAlpha',0.7);
h3 = histogram(all_spots,'BinWidth', .1, 'FaceColor','none','EdgeColor','k', 'LineWidth', 4, 'DisplayStyle', 'stairs')

% set(gca,'yscale','log')
% figure;
% histfit(all_spots, 40, 'normal')
% pd = fitdist(all_spots','Kernel','Kernel','normal', 'Support', 'positive', 'Width', .015);
pd = fitgmdist(double(all_spots'),2);
hold on
x_values = min(all_spots):.01:max(all_spots);
y = pdf(pd,x_values');
% plot(x_values,y, 'LineWidth', 3)
% set(gca,'yscale','log')
% xlim([4.07, 6.7])
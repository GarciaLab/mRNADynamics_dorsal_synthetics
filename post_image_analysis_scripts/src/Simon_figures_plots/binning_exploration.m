function binning_exploration

%% load stuff

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat']);%, 'dorsalResultsDatabase')
AllNC12Struct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);

Datasets = {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3','1DgVW','TwiPE','1Dg-8D','1Dg-5','1Dg-12'};


%% FIRST EXPLORATION moving averaging window.
% first we sort all nuclei by Dorsal concentration. Then we calculate a
% moving average of Dorsal concentration across a group of nuclei and a
% moving average of some transcriptional metric for that group. We move
% this window one nucleus at a time and plot

%IMPORTANT NOTE: we are combining 1X (FFF) with 2X.

clearvars -except AllNC12Struct Datasets

AllNC12Table = struct2table(AllNC12Struct); % convert the struct to a table
sortedT = sortrows(AllNC12Table, 'dorsalFluoFeature'); % sort the table by dorsal fluo
sortedStruct = table2struct(sortedT); % change it back to struct array 
pseudoBinN = 40;
count=1;

figure
tiledlayout('flow')
hold on
for e = 1:length(Datasets)
    Datasets{e};
    EnhancerStruct = sortedStruct(contains({sortedStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
    DorsalFluos = [];
    OutputmRNA = [];
    count1=1;
    for i = 1:length(EnhancerStruct)
        Dorsal = EnhancerStruct(i).dorsalFluoFeature;
        mRNA = EnhancerStruct(i).particleAccumulatedFluo;
        %mRNA = EnhancerStruct(i).particleFluo95;
%         if ~isempty(mRNA)
%             OutputmRNA(count1) = mean(mRNA);
%             DorsalFluos(count1) = Dorsal;
%             count1=count1+1;
%         end
        if ~isempty(mRNA)
            OutputmRNA(count1) = nanmean(mRNA);
            %OutputmRNA(count1) = 1;

        else
            OutputmRNA(count1) = 0;
            %OutputmRNA(count1) = NaN;

        end
        DorsalFluos(count1) = Dorsal;
        count1=count1+1;
        %plot(DorsalFluos,OutputmRNA,'o')
    end
    OutputmRNA = OutputmRNA(~isnan(DorsalFluos));
    DorsalFluos = DorsalFluos(~isnan(DorsalFluos));
    
    OutputmRNA2 = [];
    DorsalFluos2 = [];
    NuclearCounts = [];
    OutputmRNA2Error=[];
    count2 = 1;
    
    
    for i = 1:length(DorsalFluos)-pseudoBinN
        DorsalFluos2(count2) = nanmean(DorsalFluos(i:i+pseudoBinN));
        %DorsalFluos2(count2) = nanmean(DorsalFluos(i:min(i+pseudoBinN,length(DorsalFluos))));

        OutputmRNA2(count2) = nanmean(OutputmRNA(i:i+pseudoBinN));
        %OutputmRNA2(count2) = nanmean(OutputmRNA(i:min(i+pseudoBinN,length(DorsalFluos))));

%         OutputmRNA2Error(count2) = nanstd(OutputmRNA(i:i+pseudoBinN));
%         OutputmRNA2Error(count2) = OutputmRNA2Error(count2)/(pseudoBinN-1);
%         NuclearCounts(count2) = length((DorsalFluos(i:min(i+pseudoBinN,length(DorsalFluos)))));       
        count2=count2+1;
    end
    %OutputmRNA2 = cumsum(OutputmRNA2);


    nexttile
    errorbar(DorsalFluos2,OutputmRNA2,OutputmRNA2Error,'CapSize',0)
    ylim([0 500])
    xlim([0 3800])
    legend(Datasets{e},'Location','northwest')

end
%legend(Datasets)
hold off
ylabel('fraction active')
xlabel('Dorsal fluorescence (AU)')



%%

% %% sort the struct according to nucleus fluorescence
% AllNC12Table = struct2table(AllNC12Struct); % convert the struct to a table
% sortedT = sortrows(AllNC12Table, 'dorsalFluoFeature'); % sort the table by dorsal fluo
% sortedStruct = table2struct(sortedT); % change it back to struct array 
% pseudoBinN = 15;
% count=1;
% figure
% hold on
% for e = 1:length(Datasets)
%     
%     EnhancerStruct = sortedStruct(contains({sortedStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
%     DorsalFluos = [];
%     OutputmRNA = [];
%     
%     for i = 1:length(EnhancerStruct)
%         Dorsal = EnhancerStruct(i).dorsalFluoFeature;
%         mRNA = EnhancerStruct(i).particleFluo;
%         if ~isempty(mRNA)
%             %mRNA = mean(mRNA);
%             mRNA = 1;
%         else
%             mRNA = 0; 
%         end
%         DorsalFluos(i) = Dorsal;
%         OutputmRNA(i) = mRNA;
%     end
%     
%     % here we correct for the number of nuclei
%     DorsalFluos2=[];
%     OutputmRNA2 = [];
%     count=1;
%     for i = 0:pseudoBinN:length(DorsalFluos)
%         DorsalFluos2(count) = nanmean(DorsalFluos(i+1:min(i+pseudoBinN,length(DorsalFluos))));
%         OutputmRNA2(count) = sum(OutputmRNA(i+1:min(i+pseudoBinN,length(DorsalFluos)))); 
%         count=count+1;
%     end
%     
% %     CumDist = cumsum(OutputmRNA);
%     CumDist2 = cumsum(OutputmRNA2);
%     %CumDist = CumDist./CumDist(end);
%     CumDist2 = CumDist2./CumDist2(end);
% 
%     plot(DorsalFluos2,CumDist2,'-','LineWidth',2)
% end
% legend(Datasets)
% hold off
% ylabel('dNactive/d[Dl]')
% xlabel('Dorsal fluorescence (AU)')
% 
% 
% %% same but only consider nuclei with spots
% AllNC12Table = struct2table(AllNC12Struct); % convert the struct to a table
% sortedT = sortrows(AllNC12Table, 'dorsalFluoFeature'); % sort the table by dorsal fluo
% sortedStruct = table2struct(sortedT); % change it back to struct array 
% pseudoBinN = 15;
% count=1;
% figure
% hold on
% for e = 1:length(Datasets)   
%     EnhancerStruct = sortedStruct(contains({sortedStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
%     DorsalFluos = [];
%     OutputmRNA = [];
%     count1=1;
%     for i = 1:length(EnhancerStruct)
%         Dorsal = EnhancerStruct(i).dorsalFluoFeature;
%         mRNA = EnhancerStruct(i).particleFluo;
%         if ~isempty(mRNA)
%             OutputmRNA(count1) = mean(mRNA);
%             DorsalFluos(count1) = Dorsal;
%             count1=count1+1;
%         end
%     end
%     
%     
%     % correct for the number of nuclei
%     DorsalFluos2=[];
%     OutputmRNA2 = [];
%     count2=1;
%     for i = 0:pseudoBinN:length(DorsalFluos)
%         DorsalFluos2(count2) = nanmean(DorsalFluos(i+1:min(i+pseudoBinN,length(DorsalFluos))));
%         OutputmRNA2(count2) = sum(OutputmRNA(i+1:min(i+pseudoBinN,length(DorsalFluos)))); 
%         count2=count2+1;
%     end
%     
% %     CumDist = cumsum(OutputmRNA);
%     CumDist2 = cumsum(OutputmRNA2);
%     %CumDist = CumDist./CumDist(end);
%     %CumDist2 = CumDist2./CumDist2(end);
% 
%     plot(DorsalFluos2,CumDist2,'-','LineWidth',2)
% end
% legend(Datasets)
% hold off
% ylabel('dNactive/d[Dl]')
% xlabel('Dorsal fluorescence (AU)')
% 
% 
% 
% %% Bin Dorsal such that there are linear increments in the accumulated mRNA
% Nbins = 15;
% FluoBinsValues = [0];
% for c = linspace(1/Nbins,1,Nbins)
%     % find position in the CDF array that is closest to this bin
%     [~,idx] = min(abs(CumDist2-c));
%     % find what Dl concentration this corresponds to
%     FluoBinsValues = [FluoBinsValues DorsalFluos2(idx)];
% end
% FluoBinsValues
% 
% figure
% plot(DorsalFluos2,CumDist2,'-ko','LineWidth',1)
% hold on
% for b = FluoBinsValues
%     plot([b b],[0 1],'r-')
% end
% hold off
% 
% ylabel('integrated mRNA (normalized to max)')
% xlabel('Dorsal fluorescence (AU)')

%% moving window of Dorsal fluorescence
clearvars -except AllNC12Struct Datasets

AllNC12Table = struct2table(AllNC12Struct); % convert the struct to a table
sortedT = sortrows(AllNC12Table, 'dorsalFluoFeature'); % sort the table by dorsal fluo
sortedStruct = table2struct(sortedT); % change it back to struct array 
pseudoBinN = 60;
count=1;

figure
tiledlayout('flow')
hold on
for e = 1:length(Datasets)
    Datasets{e}
    EnhancerStruct = sortedStruct(contains({sortedStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
    DorsalFluos = [];
    OutputmRNA = [];
    count1=1;
    for i = 1:length(EnhancerStruct)
        Dorsal = EnhancerStruct(i).dorsalFluoFeature;
        mRNA = EnhancerStruct(i).particleAccumulatedFluo;
        %mRNA = EnhancerStruct(i).particleFluo95;
%         if ~isempty(mRNA)
%             OutputmRNA(count1) = mean(mRNA);
%             DorsalFluos(count1) = Dorsal;
%             count1=count1+1;
%         end
        if ~isempty(mRNA)
            OutputmRNA(count1) = nanmean(mRNA);
            %OutputmRNA(count1) = 1;

        else
            OutputmRNA(count1) = 0;
            %OutputmRNA(count1) = NaN;

        end
        DorsalFluos(count1) = Dorsal;
        count1=count1+1;
        %plot(DorsalFluos,OutputmRNA,'o')
    end
    OutputmRNA = OutputmRNA(~isnan(DorsalFluos));
    DorsalFluos = DorsalFluos(~isnan(DorsalFluos));
    
    OutputmRNA2 = [];
    DorsalFluos2 = [];
    NuclearCounts = [];
    OutputmRNA2Error=[];
    count2 = 1;
    
    
    for i = 1:length(DorsalFluos)-pseudoBinN
        DorsalFluos2(count2) = nanmean(DorsalFluos(i:i+pseudoBinN));
        %DorsalFluos2(count2) = nanmean(DorsalFluos(i:min(i+pseudoBinN,length(DorsalFluos))));

        OutputmRNA2(count2) = nanmean(OutputmRNA(i:i+pseudoBinN));
        %OutputmRNA2(count2) = nanmean(OutputmRNA(i:min(i+pseudoBinN,length(DorsalFluos))));

%         OutputmRNA2Error(count2) = nanstd(OutputmRNA(i:i+pseudoBinN));
%         OutputmRNA2Error(count2) = OutputmRNA2Error(count2)/(pseudoBinN-1);
%         NuclearCounts(count2) = length((DorsalFluos(i:min(i+pseudoBinN,length(DorsalFluos)))));       
        count2=count2+1;
    end
    %OutputmRNA2 = cumsum(OutputmRNA2);


    nexttile
    errorbar(DorsalFluos2,OutputmRNA2,OutputmRNA2Error,'CapSize',0)
    %ylim([0 500])
    xlim([0 3800])
    ylim([0 400])
    legend(Datasets{e},'Location','northwest')

end
%legend(Datasets)
hold off
ylabel('accumulated fluorescence')
xlabel('Dorsal fluorescence (AU)')



%% show every single nucleus
clearvars -except AllNC12Struct Datasets

AllNC12Table = struct2table(AllNC12Struct); % convert the struct to a table
sortedT = sortrows(AllNC12Table, 'dorsalFluoFeature'); % sort the table by dorsal fluo
sortedStruct = table2struct(sortedT); % change it back to struct array 
%pseudoBinN = 40;
count=1;
figure
tiledlayout('flow')
hold on
for e = 1:length(Datasets)
    Datasets{e}
    EnhancerStruct = sortedStruct(contains({sortedStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
    DorsalFluos = [];
    OutputmRNA = [];
    count1=1;
    for i = 1:length(EnhancerStruct)
        Dorsal = EnhancerStruct(i).dorsalFluoFeature;
        mRNA = EnhancerStruct(i).particleAccumulatedFluo;
        %mRNA = EnhancerStruct(i).particleFluo95;
        if ~isempty(mRNA)
            %OutputmRNA(count1) = nanmean(mRNA);
            OutputmRNA(count1) = 1;
        else
            OutputmRNA(count1) = 0;
            %OutputmRNA(count1) = NaN;

        end
        DorsalFluos(count1) = Dorsal;
        count1=count1+1;

    end

    nexttile
    plot(DorsalFluos,OutputmRNA,'ko')
    ylim([-.3 1])
    xlim([0 3800])
    legend(Datasets{e})

end
%legend(Datasets)
hold off
ylabel('fraction active')
xlabel('Dorsal fluorescence (AU)')




end

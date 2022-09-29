function plotDorsalDosage

%% get the prefix names
%this is where the data is
load('S:\Simon\Dropbox\DorsalSyntheticsDropbox\dorsalResultsDatabase.mat')
condIdx = contains({combinedCompiledProjects_allEnhancers.prefix},'1Dg11_2xDl') & [combinedCompiledProjects_allEnhancers.cycle]==12;

combinedCompiledProjects_SomeEnhancers = combinedCompiledProjects_allEnhancers(condIdx);

Prefixes = string({combinedCompiledProjects_SomeEnhancers.prefix});
UniquePrefixes = unique(Prefixes);

Idx1X = zeros(1,length(UniquePrefixes));
Idx2X = zeros(1,length(UniquePrefixes));

% calculate the offset from embryos without dl-Venus
OffsetPrefixes = {'2020-10-23-dl1-MCPmCh_VenusBackground_1','2020-10-23-dl1-MCPmCh_VenusBackground_2'};
%VenusOffset = getBackgroundFluo(OffsetPrefixes);
VenusOffset = 15; %AUs

%% loop over prefixes to see what's the fluorescence of their nuclei
MeanFluoPerPrefix = [];
for p = 1:length(UniquePrefixes)
    PrefixFluoValues = [];
    Prefix = UniquePrefixes{p};
    
    % to keep track of whether it's a 1x or 2x embryo
    if contains(Prefix,'EfEfEf') || contains(Prefix,'FFF')
        Idx1X(p) = 1;
    elseif contains(Prefix,'2x') || contains(Prefix,'2X')
        Idx2X(p) = 1;
    end
    
    
    %loop over all nuclei (from all projects)
    for n = 1:length(combinedCompiledProjects_allEnhancers)
        nucleusPrefix = combinedCompiledProjects_allEnhancers(n).prefix;
        nucleusNC = combinedCompiledProjects_allEnhancers(n).cycle;
        % if the nucleus belongs to the current prefix and it's in NC12,
        % add its fluo to this prefix
        if strcmpi(Prefix,nucleusPrefix) && nucleusNC == 12
            PrefixFluoValues = [PrefixFluoValues combinedCompiledProjects_allEnhancers(n).dorsalFluoFeature];
        end
    end
    % store the mean fluo of all nuclei belonging to each prefix
    MeanFluoPerPrefix = [MeanFluoPerPrefix (nanmean(PrefixFluoValues)-VenusOffset)];
end

MeanFluoPer1XPrefix = MeanFluoPerPrefix;
%MeanFluoPer1XPrefix(logical(Idx2X)) = nan;
MeanFluoPer2XPrefix = MeanFluoPerPrefix;
%MeanFluoPer2XPrefix(logical(Idx1X)) = nan;


%% find the top brightest (or dimmest) prefixes
TopNumber = 10; %number of prefixes to plot

sortedFluos1XPrefixes = sort(MeanFluoPer1XPrefix(find(~isnan(MeanFluoPer1XPrefix))),'descend'); %default is 'ascend'
Extreme1Xidx = find(ismember(MeanFluoPer1XPrefix, sortedFluos1XPrefixes(1:TopNumber)));
Extreme1XPrefixes = UniquePrefixes(Extreme1Xidx);

sortedFluos2XPrefixes = sort(MeanFluoPer2XPrefix(find(~isnan(MeanFluoPer2XPrefix))),'descend'); %default is ascending
Extreme2Xidx = find(ismember(MeanFluoPer2XPrefix, sortedFluos2XPrefixes(1:TopNumber)));
Extreme2XPrefixes = UniquePrefixes(Extreme2Xidx);


%%
TimeSinceAnaphase = [0:0.4:12]; % minutes
AllNuclearTraces1X = nan(1,length(TimeSinceAnaphase));
MeanAllNuclearTraces1X = nan(TopNumber,length(TimeSinceAnaphase));
RowCounterB=1;
figure
colors = {'k','b','g','y','r','m','c',[.8 0.2 .8],[.5 .5 .5],[0.2 .8 .8],[.2 .8 .2],[.9 0.1 .1],[.1 0.9 .1]};
hold on
for p = Extreme1XPrefixes(1,[2,6,9,10])
    RowCounterB
    %ThisPrefixName = Brightest1XPrefixes{p};
    ThisPrefixNuclearTraces = nan(40,length(TimeSinceAnaphase));
    RowCounterA = 1;
    for n = 1:length(combinedCompiledProjects_allEnhancers)
        nucleusPrefix = combinedCompiledProjects_allEnhancers(n).prefix; 
        nucleusNC = combinedCompiledProjects_allEnhancers(n).cycle;
        nucleusFluo = combinedCompiledProjects_allEnhancers(n).dorsalFluoTimeTrace;
        if strcmpi(p,nucleusPrefix) && nucleusNC==12 && ~isempty(nucleusFluo)
           FluoTimeTrace = combinedCompiledProjects_allEnhancers(n).dorsalFluoTimeTrace - VenusOffset;
           AbsTimeTrace = combinedCompiledProjects_allEnhancers(n).nuclearTimeSinceAnaphase;
           plot(AbsTimeTrace,FluoTimeTrace,'Color',colors{RowCounterB});
           for t = 1:length(AbsTimeTrace)
               time = AbsTimeTrace(t);
               fluo = FluoTimeTrace(t);
               [~,column] = min(abs(TimeSinceAnaphase-time));
               ThisPrefixNuclearTraces(RowCounterA,column)= fluo;
           end
           RowCounterA = RowCounterA+1;
        end 
        MeanAllNuclearTraces1X(RowCounterB,:) = smooth(nanmean(ThisPrefixNuclearTraces));
    end
    AllNuclearTraces1X = [AllNuclearTraces1X;ThisPrefixNuclearTraces];
    plot(TimeSinceAnaphase,smooth(nanmean(ThisPrefixNuclearTraces)),'Color',colors{RowCounterB})
    RowCounterB = RowCounterB+1;
end
hold off
title('1x')
legend('1','2','3','4','5','6','7','8','9','10')


%
AllNuclearTraces2X = nan(1,length(TimeSinceAnaphase));
RowCounterB=1;
MeanAllNuclearTraces2X = nan(TopNumber,length(TimeSinceAnaphase));
figure
%colors = {'k','b','g','y','r','m','c',[.8 0.2 .8],[.5 .5 .5],[0.2 .8 .8]};
hold on
for p = Extreme2XPrefixes(1,[2,3,7,8])
    RowCounterB
    ThisPrefixNuclearTraces = nan(40,length(TimeSinceAnaphase));
    RowCounterA = 1;
    for n = 1:length(combinedCompiledProjects_allEnhancers)
        nucleusPrefix = combinedCompiledProjects_allEnhancers(n).prefix;
        nucleusNC = combinedCompiledProjects_allEnhancers(n).cycle;
        if strcmpi(p,nucleusPrefix) && nucleusNC == 12
           FluoTimeTrace = combinedCompiledProjects_allEnhancers(n).dorsalFluoTimeTrace - VenusOffset;
           AbsTimeTrace = combinedCompiledProjects_allEnhancers(n).nuclearTimeSinceAnaphase;
           %plot(AbsTimeTrace,FluoTimeTrace,'Color',colors{RowCounterB});
           %waitforbuttonpress
           for t = 1:length(AbsTimeTrace)
               time = AbsTimeTrace(t);
               fluo = FluoTimeTrace(t);
               [~,column] = min(abs(TimeSinceAnaphase-time));
               ThisPrefixNuclearTraces(RowCounterA,column)= fluo;
           end
           RowCounterA = RowCounterA+1;
        end
        MeanAllNuclearTraces2X(RowCounterB,:) = smooth(nanmean(ThisPrefixNuclearTraces));
    end
    AllNuclearTraces2X = [AllNuclearTraces2X;ThisPrefixNuclearTraces];
    plot(TimeSinceAnaphase,smooth(nanmean(ThisPrefixNuclearTraces)),'Color',colors{RowCounterB})
    RowCounterB = RowCounterB+1;
end
hold off
title('2x')
legend('1','2','3','4','5','6','7','8','9','10')

AllNuclearTraces1X(AllNuclearTraces1X==0)=nan;
AllNuclearTraces2X(AllNuclearTraces2X==0)=nan;
NnucleiPerTime1X = sum(~isnan(AllNuclearTraces1X));
NnucleiPerTime2X = sum(~isnan(AllNuclearTraces2X));

%% Plots

figure
moleculesPerAU = 1/6; %this is from getVenusCalibrationSA([],[])
backgroundFluo = 50; %this was measured in embryos lacking Dorsal-Venus but having mCherry
yyaxis left
hold on
errorbar(TimeSinceAnaphase,nanmean(AllNuclearTraces1X)-backgroundFluo,nanstd(AllNuclearTraces1X)./NnucleiPerTime1X,'-r')
errorbar(TimeSinceAnaphase,nanmean(AllNuclearTraces2X)-backgroundFluo,nanstd(AllNuclearTraces2X)./NnucleiPerTime1X,'-b')
ylabel('nuclear Dorsal fluorescence(AU)')
ylim([0 3500])
yyaxis right
plot(TimeSinceAnaphase,(nanmean(AllNuclearTraces1X)-backgroundFluo).*(moleculesPerAU/2),'LineStyle','none')
ylabel('Dorsal dimers concentration (nM)')
ylim([0 3500].*(moleculesPerAU/2))
xlim([0 12])
legend('1X','2X')


figure
dosageFactor = 1;%1.5;
yyaxis left
hold on
errorbar(TimeSinceAnaphase,nanmean(MeanAllNuclearTraces1X)-backgroundFluo,nanstd(MeanAllNuclearTraces1X)./4,...
    '-ro','CapSize',0,'MarkerFaceColor','r','MarkerEdgeColor','none')
errorbar(TimeSinceAnaphase,nanmean(MeanAllNuclearTraces2X)-backgroundFluo,nanstd(MeanAllNuclearTraces2X)./4,...
    '-bo','CapSize',0,'MarkerFaceColor','b','MarkerEdgeColor','none')
plot(TimeSinceAnaphase,(nanmean(MeanAllNuclearTraces1X)-backgroundFluo).*dosageFactor,'k-')
ylabel('nuclear Dorsal fluorescence(AU)')
ylim([0 3500])
hold off
yyaxis right
plot(TimeSinceAnaphase,nanmean(MeanAllNuclearTraces1X).*(moleculesPerAU/2),'LineStyle','none')
ylabel('Dorsal dimers concentration (nM)')
ylim([0 3500].*(moleculesPerAU/2))
xlim([0 12])
legend('1X','2X')

figure
yyaxis left
hold on
shadedErrorBar(TimeSinceAnaphase,nanmean(MeanAllNuclearTraces1X)-backgroundFluo,...
    nanstd(MeanAllNuclearTraces1X)./4,'lineProps',{'Color','r','LineWidth',2})
shadedErrorBar(TimeSinceAnaphase,nanmean(MeanAllNuclearTraces2X)-backgroundFluo,...
    nanstd(MeanAllNuclearTraces2X)./4,'lineProps',{'Color','b','LineWidth',2})
ylabel('nuclear Dorsal fluorescence(AU)')
ylim([0 3500])
hold off
yyaxis right
plot(TimeSinceAnaphase,nanmean(MeanAllNuclearTraces1X).*(moleculesPerAU/2),'LineStyle','none')
ylabel('Dorsal dimers concentration (nM)')
ylim([0 3500].*(moleculesPerAU/2))
xlim([0 12])
legend('1X','2X')
 
    
end

% 
% 
% %%
% 
% 
% 
% 
% 
% 
% mothers = string([dorsalResultsDatabase.mother]);
% prefixes1X = mothers=="FFF";
% prefixes2X = mothers=="2x";
% 
% 
% 
% 
% %%
% 
% 
% 
% 
% 
% 
% % FluoFeature1_1X = [];
% % FluoFeature2_1X = [];
% % FluoFeature1_2X= [];
% % FluoFeature2_2X= [];
% Idx1X = [];
% Idx2X = [];
% DosageIdx = zeros(1,length(combinedCompiledProjects_allEnhancers));
% 
% for nucleus = 1:length(combinedCompiledProjects_allEnhancers)
%     nucleus/length(combinedCompiledProjects_allEnhancers)
%     
%     prefix = combinedCompiledProjects_allEnhancers(nucleus).prefix;
% %     fluoParameter1 = combinedCompiledProjects_allEnhancers(nucleus).dorsalFluoFeature;
% %     fluoParameter2 = nanmean(combinedCompiledProjects_allEnhancers(nucleus).dorsalFluoTimeTrace);
% %     
%     if contains(prefix,'EfEfEf')
%         Idx1X = [Idx1X nucleus];
%         DosageIdx(nucleus)=1;
%         
%     elseif contains(prefix,'2x') || contains(prefix,'2X')
%         Idx2X = [Idx2X nucleus];
%         DosageIdx(nucleus)=2;
%     end
% 
% end
% 
% %%
% AllFluos = [combinedCompiledProjects_allEnhancers.dorsalFluoFeature];
% 
% Fluo1X = AllFluos(Idx1X);
% Fluo2X = AllFluos(Idx2X);
% 
% figure
% hold on
% histogram(Fluo1X,'normalization','probability')
% histogram(Fluo2X,'normalization','probability')
% hold off
% 
% %%
% cutoff = 95;
% Top_1X = prctile(Fluo1X,cutoff);
% Top_2X = prctile(Fluo2X,cutoff);
% 
% BrighterThanTop1X = AllFluos>Top_1X;
% BrighterThanTop2X = AllFluos>Top_2X;
% 
% BrighterThanTop1X_1Xs = BrighterThanTop1X(Idx1X);
% BrighterThanTop2X_2Xs = BrighterThanTop2X(Idx2X);
% 
% Brigthest1XIdx = Idx1X(BrighterThanTop1X_1Xs);
% Brigthest2XIdx = Idx2X(BrighterThanTop2X_2Xs);
% 
% Brightest1XFluos = AllFluos(Brigthest1XIdx);
% Brightest1XFluos = Brightest1XFluos(Brightest1XFluos<4000);
% Brightest2XFluos = AllFluos(Brigthest2XIdx);
% 
% figure
% hold on
% histogram(Brightest1XFluos,'normalization','probability')
% histogram(Brightest2XFluos,'normalization','probability')
% hold off
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 
% % figure
% % hold on
% % histogram(FluoFeature1_1X,'normalization','probability')
% % histogram(FluoFeature1_2X,'normalization','probability')
% % hold off
% 
% % figure
% % hold on
% % histogram(FluoFeature2_1X,'normalization','probability')
% % histogram(FluoFeature2_2X,'normalization','probability')
% % hold off
% 
% 
% %%
% % plot fluo over time of the top x percent
% 
% %remove weird outliers from 1x data
% FluoFeature1_1X = FluoFeature1_1X(FluoFeature1_1X<3800);
% 
% Top5_1X = prctile(FluoFeature1_1X,90);
% Top5_2X = prctile(FluoFeature1_2X,90);
% 
% FluoFeature3_1X = [];
% FluoFeature3_2X= [];
% 
% for nucleus = 1:length(combinedCompiledProjects_allEnhancers)
%     nucleus/length(combinedCompiledProjects_allEnhancers)
%     
%     prefix = combinedCompiledProjects_allEnhancers(nucleus).prefix;
%     fluoParameter3 = combinedCompiledProjects_allEnhancers(nucleus).dorsalFluoFeature;
%     %fluoParameter2 = nanmean(combinedCompiledProjects_allEnhancers(nucleus).dorsalFluoTimeTrace);
%     
%     if contains(prefix,'EfEfEf') && (fluoParameter3  > Top5_1X)
%         FluoFeature3_1X = [FluoFeature3_1X  fluoParameter3];
%         %FluoFeature2_1X = [FluoFeature2_1X  fluoParameter2];
%         
%     elseif contains(prefix,'2x') || contains(prefix,'2X') && (fluoParameter3 > Top5_2X)
%         FluoFeature3_2X = [FluoFeature3_2X  fluoParameter3];
%         %FluoFeature2_2X = [FluoFeature2_2X  fluoParameter2];
%     end
% 
% end
% 
% figure
% hold on
% histogram(FluoFeature3_1X,'normalization','probability')
% histogram(FluoFeature3_2X,'normalization','probability')
% hold off
% 
% 
% 

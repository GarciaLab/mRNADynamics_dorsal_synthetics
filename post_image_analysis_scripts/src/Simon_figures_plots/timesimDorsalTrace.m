function timesimDorsalTrace(numBins,control)

DorsalFluoStruct = getDorsalFluoTraces(numBins,[]);
time_vec = DorsalFluoStruct(1).absoluteTime/60;
binValues = DorsalFluoStruct(1).binValues;
%close all

% for bin = 1:length(DorsalFluoStruct)
%    
%     dls = DorsalFluoStruct(bin).meanDorsalFluo;
% 
%     [fraction_actives(bin), mean_onsets(bin)] = timesim(time_vec,dls,[]);
% 
% end
% 
% figure
% yyaxis left
% plot(linspace(0,4500,length(fraction_actives)),fraction_actives,'b','LineWidth',2, 'DisplayName', 'fraction')
% ylabel('fraction active')
% ylim([0 1])
% 
% yyaxis right
% plot(linspace(0,4500,length(fraction_actives)),mean_onsets,'r','LineWidth',2, 'DisplayName', 'onset')
% ylabel('onset time')
% ylim([0 10])
% 
% legend()

%% now loop over Kds
% KDs = [15000,10000,8500,7000,5000,3000,2500,1000,500,100];
KDs = 10:1000:10000;
fraction_actives = nan(length(KDs),length(DorsalFluoStruct));
mean_onsets = nan(length(KDs),length(DorsalFluoStruct));

%for rep = 1:100
    counter=1;
    for k = KDs
        for bin = 1:length(DorsalFluoStruct)
            
            if isempty(control)
                dls = DorsalFluoStruct(bin).meanDorsalFluo;
            elseif control==1
                dls = ones(length(DorsalFluoStruct(bin).meanDorsalFluo),1).*DorsalFluoStruct(bin).originalFluoFeature;
            end
            
            [fraction_actives(counter,bin), mean_onsets(counter,bin)] = timesim_interp(time_vec,dls,'kd', k, 'c', 1.5, 'nOffStates', 5, 't_cycle', 8);
        end
        counter=counter+1;
    end
%end

% average across simulation replicates
% fraction_actives = nanmean(fraction_actives,3);
% mean_onsets = nanmean(mean_onsets,3);

% plots
Palette = viridis(length(KDs));
XvalsForPlot = linspace(0,4500,length(fraction_actives));

figure
hold on
for k = 1:length(KDs)
    plot(XvalsForPlot,fraction_actives(k,:),'Color',Palette(k,:),'LineWidth',2, 'DisplayName', 'fraction')
end
ylabel('fraction active')
ylim([0 1])

figure
hold on
for k = 1:length(KDs)
    plot(XvalsForPlot,mean_onsets(k,:),'Color',Palette(k,:),'LineWidth',2, 'DisplayName', num2str(KDs(k)))
end
ylabel('onset')
ylim([0 9]);
xlim([100, 4000]);
legend()

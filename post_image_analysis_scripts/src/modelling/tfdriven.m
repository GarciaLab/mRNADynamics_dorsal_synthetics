function tfdriven()
close all;

nPlots = 10;
dmax = 5000;
d = logspace(1, log10(dmax)); %aus
t = linspace(0, 10)'; %mins
kd = 500;
kds = logspace(2, 4, nPlots);
cs = logspace(-1, 2, nPlots);
R = 500;
nSteps = 5;
t_nc13 = 10;
cmap = colormap(viridis(nPlots));


%%
%functions
tonfun  = @(d, c, kd) 5.*c.^(-1).*(1+d.^(-1).*kd+2500.*c.^5.*d.^4.*(3.*d.^4+30.*c.*d.^4+ ...
    150.*c.^2.*d.^4+500.*c.^3.*d.^4+1250.*c.^4.*d.^4+12.*d.^3.*kd+90.* ...
    c.*d.^3.*kd+300.*c.^2.*d.^3.*kd+500.*c.^3.*d.^3.*kd+18.*d.^2.* ...
    kd.^2+90.*c.*d.^2.*kd.^2+150.*c.^2.*d.^2.*kd.^2+12.*d.*kd.^3+30.* ...
    c.*d.*kd.^3+3.*kd.^4+(-3).*exp(1).^(10.*c.*d.*(d+kd).^(-1)).*(d+ ...
    kd).^4).^(-1));
%% setup parameters and stuff
close all;

nPlots = 100;
cmap = colormap(viridis(nPlots));
% nPlots = 10;
d = logspace(1, 4); % range of dorsal concentrations in AUs
t = linspace(0, 10)'; % time in mins
kd = 500; % Kd of Dorsal binding in AUs
cs = linspace(5, 15, nPlots);  % rate increase per unit dorsal
% cs = logspace(-1, 2, nPlots);  % rate increase per unit dorsal
R = 1; % rate of transcription when the dorsal site is 100% occupied
nSteps = 5; % number of irreversible steps
t_nc13 = 10; % duration of the transcriptional window in minutes

% [D, T] = meshgrid(d,t); %?
% 
% figure;
% t_surf = tiledlayout('flow');

f2 = figure; ax2 = axes(f2);
f3 = figure; ax3 = axes(f3);
f4 = figure; ax4 = axes(f4);
f5 = figure; ax5 = axes(f5);

% initialize arrays containing the model outputs
dmrnadt = [];
paccessible  = [];
n = 0; %for indexing
fon = [];
ton2 = [];
cdfg = [];
temp8 = [];
cdfg2 = [];
ton3 = [];
ton4 = [];

% loop over values of c
for c = cs
    
    n = n + 1; %for indexing
    
%     %five irreversible steps, SA: I think this is a wrong equation
%     temp = (-1/24).*d.*exp(-c.*d.*t).*(d+kd).^(-1).*R.*(24+(-24).* ...
%       exp(c.*d.*t)+c.*d.*t.*(24+c.*d.*t.*(12+c.*d.*t.*(4+c.*d.*t))) ...
%       );
%      dmrnadt(n, :, :) = temp;
%      mrna(n,:) = trapz(t, temp, 1);
%  
%  temp2 = -((d R (-24 t + 
%     E^(-c d t) (-(120/(c d)) - 96 t - 36 c d t^2 - 8 c^2 d^2 t^3 - 
%        c^3 d^3 t^4)))/(24 (d + kd)))
%   

%     % probability of being in the last state as a function of [Dl] and time
%     % for a given c
%     temp3 = (1/24)*exp(-c.*d.*t).* (-24 + 24*exp(c.*d.*t) - 24.*c.*d.*t - 12.*(c.*d.*t).^2 - 4*(c.*d.*t).^3 - (c.*d.*t).^4);
%     paccessible(n, :, :) = temp3;
%     
%     %odds of reaching the end of the cycle without turning on
%     fon(n, :) = squeeze(paccessible(n, t_nc13, :)); 

    temp6 = (1/24).*d.*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-1).* ...
    R.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+c.*d.*(d+kd).^(-4).* ...
    t.*((-24).*(d+kd).^3+(-12).*c.*d.*(d+kd).^2.*t+(-4).*c.^2.*d.^2.*( ...
    d+kd).*t.^2+(-1).*c.^3.*d.^3.*t.^3));

    dmrnadt(n, :, :) = temp6;
    mrna2(n,:) = trapz(t, temp6, 1);

    % probability of being in the last state as a function of [Dl] and time
    % for a given c
    temp7 = (1/24).*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-4).*(24.* ...
      ((-1)+exp(1).^(c.*d.*(d+kd).^(-1).*t)).*kd.^4+24.*d.*kd.^3.*((-4)+ ...
      4.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t)+12.*d.^2.*kd.^2.*(( ...
      -12)+12.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*(6+c.*t))+4.* ...
      d.^3.*kd.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*( ...
      18+c.*t.*(6+c.*t)))+d.^4.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).* ...
      t)+(-1).*c.*t.*(24+c.*t.*(12+c.*t.*(4+c.*t)))));

    paccessible2(n, :, :) = temp7;   
    %odds of reaching the end of the cycle without turning on
    fon2(n, :) = squeeze(paccessible2(n, t_nc13, :));

    temp10 = 5.*c.^(-1).*(1+d.^(-1).*kd+2500.*c.^5.*d.^4.*(3.*d.^4+30.*c.*d.^4+ ...
      150.*c.^2.*d.^4+500.*c.^3.*d.^4+1250.*c.^4.*d.^4+12.*d.^3.*kd+90.* ...
      c.*d.^3.*kd+300.*c.^2.*d.^3.*kd+500.*c.^3.*d.^3.*kd+18.*d.^2.* ...
      kd.^2+90.*c.*d.^2.*kd.^2+150.*c.^2.*d.^2.*kd.^2+12.*d.*kd.^3+30.* ...
      c.*d.*kd.^3+3.*kd.^4+(-3).*exp(1).^(10.*c.*d.*(d+kd).^(-1)).*(d+ ...
      kd).^4).^(-1));

     ton4(n, :) = temp10;

     
    % plot results
    
    % 3D plots
%     nexttile(t_surf)
%     surf(D, T, squeeze(dmrnadt(n, :, :)));
%     xlim([0, 1E4]);
%     title(num2str(c));
    
%     % fraction on as a function of [Dl]
%     plot(ax2, d, fon(n, :))
%     xlabel(ax2,'Dorsal concentration (AU)')
%     ylabel(ax2,'fraction of active nuclei')
%     hold(ax2, 'on')

    % accumulated mRNA
    plot(ax4, d, mrna2(n, :))
    xlabel(ax4,'Dorsal concentration (AU)')
    ylabel(ax4,'acumulated mRNA (AU)')
    hold(ax4, 'on')
    
    
    % time on
%     plot(ax8, d, ton4(n, :))
%     xlabel(ax8,'Dorsal concentration (AU)')
%     ylabel(ax8,'time on')
%     hold(ax8, 'on')

%        
   temp3 = (1/24)*exp(-c.*d.*t).* (-24 + 24*exp(c.*d.*t) - 24.*c.*d.*t - 12.*(c.*d.*t).^2 - 4*(c.*d.*t).^3 - (c.*d.*t).^4);
    
   paccessible(n, :, :) = temp3;
   
   fon(n, :) = squeeze(paccessible(n, t_nc13, :)); %odds of reaching the end of the cycle without turning on
   
 
  
   
  
   
% mrna(n,:) = trapz(t, temp, 1);
 
 temp6 = (1/24).*d.*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-1).* ...
  R.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+c.*d.*(d+kd).^(-4).* ...
  t.*((-24).*(d+kd).^3+(-12).*c.*d.*(d+kd).^2.*t+(-4).*c.^2.*d.^2.*( ...
  d+kd).*t.^2+(-1).*c.^3.*d.^3.*t.^3));
 
 dmrnadt(n, :, :) = temp6;
 mrna2(n,:) = trapz(t, temp6, 1);

 

temp7 = (1/24).*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-4).*(24.* ...
  ((-1)+exp(1).^(c.*d.*(d+kd).^(-1).*t)).*kd.^4+24.*d.*kd.^3.*((-4)+ ...
  4.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t)+12.*d.^2.*kd.^2.*(( ...
  -12)+12.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*(6+c.*t))+4.* ...
  d.^3.*kd.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*( ...
  18+c.*t.*(6+c.*t)))+d.^4.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).* ...
  t)+(-1).*c.*t.*(24+c.*t.*(12+c.*t.*(4+c.*t)))));
 

paccessible2(n, :, :) = temp7;   
fon2(n, :) = squeeze(paccessible2(n, t_nc13, :)); %odds of reaching the end of the cycle without turning on
 
%nexttile(t_surf)
% surf(D, T, squeeze(dmrnadt(n, :, :)));
% xlim([0, 1E4]);
% title(num2str(c));



plot(ax2, d, fon(n, :))
hold(ax2, 'on')


% plot(ax3, d, mrna(n, :))
% hold(ax3, 'on')


plot(ax4, d, mrna2(n, :), 'LineWidth', 2, 'Color', cmap(n, :))
hold(ax4, 'on')


%fraction active vs. dorsal, varied c
plot(ax5, d, fon2(n, :),  'LineWidth', 2, 'Color', cmap(n, :))
hold(ax5, 'on')
xlabel(ax5,'Dorsal concentration (AU)')
ylabel(ax5,'fraction of active nuclei')
end

%%


% 
% title(ax2, 'predicted fraction active')
% xlim(ax2, [0, 3500]);
% ylim(ax2,[0, 1]);
% leg2 = legend(ax2, num2str(round(cs', 2, 'significant')));
% title(leg2, 'c')
% xlabel(ax2,'[Dorsal] (au)')
% ylabel(ax2,'fraction active')




% title(ax3, 'predicted accumulated mRNA')
% xlim(ax3, [0, 3500]);
% leg3 = legend(ax3, num2str(round(cs', 2, 'significant')));
% title(leg3, 'c')
% xlabel(ax3,'[Dorsal] (au)')
% ylabel(ax3,'normalized accumulated mRNA')


title(ax4, 'predicted accumulated mRNA.  pi \alpha occupancy')
xlim(ax4, [10, dmax]);
leg4 = legend(ax4, num2str(round(cs', 2, 'significant')));
title(leg4, 'c')
xlabel(ax4,'[Dorsal] (au)')
ylabel(ax4,'accumulated mRNA (au)')
set(ax4, 'XScale', 'log')


title(ax5, 'predicted fraction active. pi \alpha occupancy')
xlim(ax5, [10, 3500]);
ylim(ax5,[0, 1]);
% leg5 = legend(ax5, num2str(round(cs', 2, 'significant')));
% title(leg5, 'c')
xlabel(ax5,'[Dorsal] (au)')
ylabel(ax5,'fraction active')
% set(ax5, 'XScale', 'log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


figure;
c = 10;

ton4 = [];
for k = 1:length(kds)
    ton4(k, :) = tonfun(d, c, kds(k));
 plot(d, ton4(k, :), 'LineWidth', 2, 'Color', cmap(k, :))
hold on;
end

% title(ax7, 'predicted \langle T_{on} \rangleT_{on} (min)')
% xlim(ax7, [0, 3500]);
% % leg7 = legend(ax7, num2str(round(cs', 2, 'significant')));
% % title(leg7, 'c')
% xlabel(ax7,'KD (au)')
% ylabel(ax7,'mean turn on time (min)')


title('predicted T_{on} (min)')
xlim([10, dmax]);
set(gca, 'XScale', 'log')
ylim([0, 10]);
leg = legend(num2str(round(kds', 2, 'significant')));
title(leg, 'K_D')
xlabel('[Dl] (au)')
ylabel('mean turn on time (min)')

% title(ax8, 'predicted T_{on} (min)')
% xlim(ax8, [0, 3500]);
% leg8 = legend(ax8, num2str(round(cs', 2, 'significant')));
% title(leg8, 'c')
% xlabel(ax8,'[Dl] (au)')
% ylabel(ax8,'mean turn on time (min)')

%% Compare with data

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])%, 'dorsalResultsDatabase')
AllStruct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);
b = load('S:\Simon\Dropbox\DorsalSyntheticsDropbox\dorsalResultsDatabase.mat');

Datasets = {'1Dg11'}%, '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3','1DgVW'};
PatserScores = [6.23,5.81,5.39,5.13,4.8,4.73,4.29];
AllTurnOnData = struct('Embryos',[]);
Nbins = 19; %number of Dorsal concentration bins
NNuclei = 25; %maximum number of nuclei in nc12

%gather the data
for e = 1:length(Datasets) %loop over enhancers
    EnhancerStruct = AllStruct(contains({AllStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
    EnhancerPrefixes = unique({EnhancerStruct.prefix}); %get the name of the prefixes of this enhancer
    EnhancerEmbryos = struct('embryoTurnOnTimes',[]); %make a struct to store data, 'embryoTurnOnTimes' contains results for a single movie
   for p = 1:length(EnhancerPrefixes)
        PrefixTurnOnTimes = nan(Nbins,NNuclei); % nan array of dorsal bins x nuclei to populate with turn on times
        % make a substruct containing only the nuclei belonging to this prefix and containing Dl bin info
        PrefixStruct = EnhancerStruct(strcmpi({EnhancerStruct.prefix},EnhancerPrefixes{p}) &...
       ~isnan([EnhancerStruct.dorsalFluoBin]));
        % get the indices of the nuclei containing a particle
        onNucleiIdx = find(arrayfun(@(PrefixStruct)~isempty(PrefixStruct.particleTimeOn),PrefixStruct));
        % make a substruct containing only on nuclei of this prefix of this enhancer
        PrefixOnNucleiStruct = PrefixStruct(onNucleiIdx);
        for n = 1:length(PrefixOnNucleiStruct) %loop over nuclei to gather Dl bin and turn on time and add to PrefixTurnOnTimes 
            PrefixTurnOnTimes(PrefixOnNucleiStruct(n).dorsalFluoBin,n) = PrefixOnNucleiStruct(n).particleTimeOn;
        end
        % add the PrefixTurnOnTimes matrix to this enhancer struct
        EnhancerEmbryos(p).embryoTurnOnTimes = PrefixTurnOnTimes;
   end
   % add this enhancer struct to the struct containing all data
   AllTurnOnData(e).Embryos= EnhancerEmbryos;   
end

% plot stuff
figure
hold on
TurnOnTimePerKd = nan(length(Datasets),Nbins);
for enh = 1:length(Datasets)   
    EnhancerEmbryos = AllTurnOnData(enh).Embryos;
    AllsingleEmbryoMeans = nan(1,19);
    for emb = 1:length(EnhancerEmbryos)
        embryoTurnOnTimes = EnhancerEmbryos(emb).embryoTurnOnTimes;
        singleEmbryoMean = nanmean(embryoTurnOnTimes,2); %mean turn on across nuclei in the same Dl bin
        AllsingleEmbryoMeans(emb,:) = singleEmbryoMean; % make an array containing the means of each embryo of this enhancer
    end
    NEmbryos = sum(~isnan(AllsingleEmbryoMeans));
    mean = nanmean(AllsingleEmbryoMeans);
    SEM = nanstd(AllsingleEmbryoMeans,1)./NEmbryos;
    mean(mean>10) = nan;
    errorbar([0:125:3750],mean,SEM,'LineWidth',2,'CapSize',0)
    TurnOnTimePerKd(enh,:) = mean;
    meanTurnOnTimePerKd(enh) = nanmean(mean);
    SEMturnOnTimePerKd(enh) = nanstd(mean)./sum(~isnan(mean));
end
hold off
legend(Datasets)
xlabel('Dorsal concentration')
ylabel('mean turn on time')
ylim([0 10])
set(gca,'XScale','log')


figure
hold on
plot(PatserScores,TurnOnTimePerKd,'o','MarkerFaceColor','r','MarkerEdgeColor','none')
errorbar(PatserScores,meanTurnOnTimePerKd,SEMturnOnTimePerKd,'ko','MarkerFaceColor','k','CapSize',0)
%ylim([0 10])
xlabel('patser score')
ylabel('mean turn on time across [Dl]')
set(gca, 'XDir','reverse')
hold off


        
   
   












end

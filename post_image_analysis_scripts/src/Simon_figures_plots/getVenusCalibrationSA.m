function getVenusCalibrationSA(MyPrefixes, BcdGFPPrefixes)


BcdGFPPrefixes = {'2019-08-15-BcdGFP-HisRFP','2019-07-26-BcdGFP-HisRFP','2019-07-08-Bcd-GFP-His-RFP',...
    '2017-07-31-BcdGFP-HisRFP','2017-07-10-BcdGFP-HisRFP','2017-06-02-BcdGFP-HisRFP','2017-06-01-BcdGFP-HisRFP',...
    '2017-05-23-BcdGFP-HisRFP',};

%% Find at what time in nc14 BcdGFP levels peak
BcdGFPPath = 'S:\Jonathan\Dropbox\ZeldaLizHernan\Data\';
NC14BcdPeakTime = []; %to store results
Palette = jet(length(BcdGFPPrefixes));
figure
for rep = 1:length(BcdGFPPrefixes)
    ReplicatePath = [BcdGFPPath BcdGFPPrefixes{rep}];
    load([ReplicatePath '\FrameInfo.mat'])
    load([ReplicatePath '\CompiledNuclei'],...
        'NParticlesAP','SDVectorAP','MeanVectorAP','APbinID','nc14');
    
    AbsTime = [FrameInfo.Time];
    TimeTrace = smooth(nanmean(MeanVectorAP'),8);
    NC14TimeTrace = TimeTrace(nc14:end);
    NC14BcdPeakFrame = find(TimeTrace == max(NC14TimeTrace));
    NC14BcdPeakTime = [NC14BcdPeakTime AbsTime(NC14BcdPeakFrame)-AbsTime(nc14)];

    % look at the results
    Color = Palette(rep,:);
    plot(AbsTime(nc14:end)-AbsTime(nc14),TimeTrace(nc14:end),'Color',Color)
    hold on
    plot([AbsTime(NC14BcdPeakFrame-nc14) AbsTime(NC14BcdPeakFrame-nc14)],[0 max(TimeTrace)],'Color',Color)
end
hold off

NC14BcdAbsPeakTime = mean(NC14BcdPeakTime); %in seconds, this is 430
NC14BcdAbsPeakTime = 430;
TimeTo15min = (15*60)-NC14BcdAbsPeakTime; %in seconds, this is 470
TimeTo15min = 470;



%% Deal with our data now
%these are our Venus-Bcd datasets
MyPrefixes = {'2020-09-05-Venus-Bcd_BcdE1_dorsalsynthetics_settings_4','2020-09-05-Venus-Bcd_BcdE1_dorsalsynthetics_settings_5',...
    '2020-09-05-Venus-Bcd_BcdE1_dorsalsynthetics_settings_7'};

% calculate the offset fluorescence in the venus channel
OffsetPrefixes = {'2020-10-23-dl1-MCPmCh_VenusBackground_1','2020-10-23-dl1-MCPmCh_VenusBackground_2'};
%VenusOffset = getBackgroundFluo(OffsetPrefixes); %this is ~15 AU
VenusOffset  =15;

DynamicsResultsPath = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox';
NReps = length(MyPrefixes);
APBins = 21;
MeanData = nan(NReps,APBins);
SDData = nan(NReps,APBins);
NData = nan(NReps,APBins);
APLengthData = nan(NReps,APBins);
%Palette = parula(NReps);
figure
Colors = {'r','g','b','k','m','c'}

for rep = 1:NReps
    ReplicatePath = [DynamicsResultsPath '/' MyPrefixes{rep}];
    load([ReplicatePath '/CompiledNuclei.mat'],...
        'NParticlesAP','SDVectorAP','MeanVectorAP','APbinID','nc14');
    load([ReplicatePath '/FrameInfo.mat'])
    AbsTime = [FrameInfo.Time];
    
    %*^*^*^*^* background subtraction and dosage correaction done here *^*
    TimeTrace = 1.5.*(smooth(nanmean(MeanVectorAP'),8) - VenusOffset);
    NC14TimeTrace = TimeTrace(nc14:end);
    NC14BcdPeakFrame = find(TimeTrace == max(NC14TimeTrace));
    
    %we use the peak Bcd fluo to calibrate absolute time into nc14
    time15minIntoNC14 = AbsTime(NC14BcdPeakFrame) + TimeTo15min;
    % the fiducial frame is 15 min into nc14.
    [~,FiducialFrame] = min(abs(AbsTime-time15minIntoNC14));
    
    %Color = Palette(rep,:);
    plot(AbsTime(nc14:end)-AbsTime(nc14),TimeTrace(nc14:end),'Color',Colors{rep})
    hold on
    plot([AbsTime(NC14BcdPeakFrame)-AbsTime(nc14) AbsTime(NC14BcdPeakFrame)-AbsTime(nc14)],[0 max(NC14TimeTrace)],'Color',Colors{rep})
    plot([AbsTime(FiducialFrame)-AbsTime(nc14) AbsTime(FiducialFrame)-AbsTime(nc14)],[0 max(NC14TimeTrace)],'Color',Colors{rep})
    ylim([0 400])
    
    APBinsCovered = ~isnan(nansum(MeanVectorAP));
    
    DosageFactor = 1;
    MeanData(rep,APBinsCovered)= DosageFactor.*(MeanVectorAP(FiducialFrame,APBinsCovered)- VenusOffset);
    SDData(rep,APBinsCovered) = sqrt(DosageFactor).*(SDVectorAP(FiducialFrame,APBinsCovered));
    NData(rep,APBinsCovered) = NParticlesAP(FiducialFrame,APBinsCovered);
    APLengthData(rep,APBinsCovered) = APbinID(APBinsCovered);
end

MeanData(MeanData==0)=nan;
hold off

figure
errorbar(APLengthData(:),MeanData(:),SDData(:),'o','CapSize',0)

% load absolute Bcd concentration data from Gregor 2007
bcd_abs_path = 'S:\Armando\Dropbox\DorsalSyntheticsDropbox\Data_Venus_calibration\';
bkg_data = readtable([bcd_abs_path 'Gregor2007Black.csv']);
bkg_data_clean.AP = bkg_data.Var1;
bkg_data_clean.nM = bkg_data.Var2;
bcd_data01 = readtable([bcd_abs_path 'Gregor2007BcdRed.csv']);
bcd_data02 = readtable([bcd_abs_path 'Gregor2007BcdBlue.csv']);


bcd_data_AP1 = [bcd_data01.Var1];
bcd_data_nM1 = [bcd_data01.Var2];
bcd_data_AP2 = [bcd_data02.Var1];
bcd_data_nM2 = [ bcd_data02.Var2];

figure
plot(bcd_data_AP1,bcd_data_nM1,'ro')
hold on
plot(bcd_data_AP2,bcd_data_nM2,'bo')
hold off

figure
yyaxis left
errorbar(APLengthData(:),MeanData(:),SDData(:),'o','CapSize',0,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8)
ylim([0 420])
ylabel('Venus-Bcd fluorescence (AU)')

yyaxis right
hold on
plot(bcd_data_AP1,bcd_data_nM1,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8)
plot(bcd_data_AP2,bcd_data_nM2,'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',8)
ylabel('EGFP-Bcd concentration (nM)')
hold off
xlabel('AP position (x/L)')



MeanDataVector = MeanData(~isnan(MeanData));
SDVector = SDData(~isnan(MeanData));%./NData(~isnan(MeanData));
APPosVector = APLengthData(~isnan(MeanData));
nM1Vector = nan(size(APPosVector));
nM2Vector = nan(size(APPosVector));

for i = 1:length(APPosVector)
    APPos = APPosVector(i);
    [~,closestGregorPos1] = min(abs(bcd_data_AP1-APPos));
    [~,closestGregorPos2] = min(abs(bcd_data_AP2-APPos));
    nM1Vector(i) = bcd_data_nM1(closestGregorPos1);
    nM2Vector(i) = bcd_data_nM2(closestGregorPos2);
end


% %% FINAL FIGURE V1, FREE FIT ALLOWING FOR Y-INTERCEPT
% figure
% hold on
% MeanDataVector(MeanDataVector==0)=nan;
% DosageCorrection = 3/2; %to account for the fact that Venus-Bcd is less abundant than EGFP-Bcd (based on cephalic furrow position)
% MeanDataVector = MeanDataVector.*DosageCorrection;
% SDVector(SDVector==0)=nan;
% errorbar(nM1Vector,MeanDataVector,SDVector,'ro','CapSize',0)
% errorbar(nM2Vector,MeanDataVector,SDVector,'bo','CapSize',0)
% xlabel('nM')
% ylabel('Venus fluorescence')
% 
% 
% % now do the fitting allowing a y intercept
% XdataVector = [nM1Vector;nM2Vector]; %concentration
% YdataVector = [MeanDataVector;MeanDataVector]; %fluorescence
% %FluoError = [SDVector;SDVector];
% 
% regX = [ones(size(XdataVector)) XdataVector];% the padding with ones is needed for the thing to run
% [b,bint,r,rint,stats] = regress(YdataVector,regX);
% slope = b(2);
% minslope = bint(2,2) ; % top of 95% confidence interval
% maxslope = bint (2,1);% bottom of 95% confidence interval
% Yintersect = b(1);
% Rsquared = stats(1);
% % add a zero to the X data for plotting purposes
% XdataVector = [0;XdataVector];
% 
% plot(XdataVector,Yintersect+XdataVector*slope,'k-','LineWidth',1)
% hold off
% ylim([0 200])
% xlim([0 45])
% title(['Slope = ' num2str(slope) '\pm' num2str(minslope-slope) ' AU/nM ; R^2 = ' num2str(Rsquared)])


%% FINAL FIGURE V2, GOING THROUGH THE ORIGIN
APPosVector = APPosVector(~isnan(MeanDataVector));
SDVector = SDVector(~isnan(MeanDataVector));
MeanDataVector = MeanDataVector(~isnan(MeanDataVector));

nM1Vector=[];nM2Vector=[];
for i = 1:length(APPosVector)
    APPos = APPosVector(i)
    [dummy,closestGregorPos1] = min(abs(bcd_data_AP1-APPos));
    [dummy,closestGregorPos2] = min(abs(bcd_data_AP2-APPos));
    nM1Vector(i) = bcd_data_nM1(closestGregorPos1);
    nM2Vector(i) = bcd_data_nM2(closestGregorPos2);
end




figure
hold on
errorbar(nM1Vector,MeanDataVector,SDVector,'ro','CapSize',0)
errorbar(nM2Vector,MeanDataVector,SDVector,'bo','CapSize',0)
xlabel('nM')
ylabel('Venus fluorescence')


% now do the fitting forcing it to go through zero
XdataVector = [nM1Vector nM2Vector]; %concentration
YdataVector = [MeanDataVector;MeanDataVector]; %fluorescence

%FluoError = [SDVector;SDVector];

slope = XdataVector'\YdataVector; % linear regression throug origin
yfit = XdataVector*slope; % calculate fitted line

% calculate R2          
f = yfit;
Bbar = mean(XdataVector');
SStot = sum((XdataVector' - Bbar).^2);
SSreg = sum((f - Bbar).^2);
SSres = sum((XdataVector' - f).^2);
R2 = 1 - SSres/SStot;
R = corrcoef(YdataVector,XdataVector');
Rsq = R(1,2).^2;

% add a zero to data for plotting purposes
XdataVector = [0 XdataVector];
yfit = [0 yfit];

plot(XdataVector,yfit,'k-','LineWidth',1)
hold off
ylim([0 270])
xlim([0 45])
%title(['Slope = ' num2str(slope) '\pm' num2str(minslope-slope) ' AU/nM ; R^2 = ' num2str(Rsquared)])
title(['Slope = ' num2str(slope) ' AU/nM ; R^2 = ' num2str(Rsq)])
   
       
    
end
    
    
% get CompiledNuclei 

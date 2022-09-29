%%
EGFPBcdData = [536.699 165.660 473.240 151.515 512.610 214.154 519.765 170.932 516.866...
    164.648 502.719 155.007 505.036 137.654 522.713 173.999 492.687 143.017];
VenusBcdData = [485.103 209.740 492.362 195.152 503.665 176.181 519.822 205.608 492.615...
    214.140 519.664 201.229 517.925 212.384 517.165 223.465 510.158 207.132 519.401 212.398];

EGFPBcdEmbryoLength = EGFPBcdData(1:2:end);
EGFPBcdCFLength = EGFPBcdData(2:2:end);
CFDataEGFP = EGFPBcdCFLength./EGFPBcdEmbryoLength;

VenusBcdEmbryoLength = VenusBcdData(1:2:end);
VenusBcdCFLength = VenusBcdData(2:2:end);
CFDataVenus = VenusBcdCFLength./VenusBcdEmbryoLength;

%% plot
figure(1)
hold on
plot(1,CFDataEGFP,'go','MarkerFaceColor','g','MarkerEdgeColor','none','MarkerSize',10)
errorbar(1,mean(CFDataEGFP),std(CFDataEGFP),'ko','CapSize',0)
plot(2,CFDataVenus,'yo','MarkerFaceColor','y','MarkerEdgeColor','none','MarkerSize',10)
errorbar(2,mean(CFDataVenus),std(CFDataVenus),'ko','CapSize',0)
plot([0.5 2.5],[0.343 0.343],'k--')
xlim([0.5 2.5])
xticks([1,2])
xticklabels({'EGFP-Bcd','Venus-Bcd'})
ylabel({'cephalic furrow position', '(% embryo length)'})

figure(2)
CFDataEGFPNorm = CFDataEGFP./0.343;
CFDataVenusNorm = CFDataVenus./0.343;
hold on
plot(1,CFDataEGFPNorm,'go','MarkerFaceColor','g','MarkerEdgeColor','none','MarkerSize',10)
errorbar(1,mean(CFDataEGFPNorm),std(CFDataEGFPNorm),'ko','CapSize',0)
plot(2,CFDataVenusNorm,'yo','MarkerFaceColor','y','MarkerEdgeColor','none','MarkerSize',10)
errorbar(2,mean(CFDataVenusNorm),std(CFDataVenusNorm),'ko','CapSize',0)
plot([0.5 2.5],[1 1],'k--')
%xlim([0.5 2.5])
xticks([1,2])
xticklabels({'EGFP-Bcd','Venus-Bcd'})
ylabel({'cephalic furrow position normalized to wt', '(% embryo length)'})

% According to Liu, the slope of the CF vs protein dosage function is 44%
% (normalized to wt)
% so an increase in the protein dosage of 100% (i.e 2x wt) leads to moving
% the cephalic furrow to 1.44 the wt.
% Conversely, a fly with a CF position X% of wt has a protein dosage of
% X/0.44

% calculate the dosage with respect to "wt"
EGFPDosage = 1-((1-mean(CFDataEGFPNorm))./0.44);
VenusDosage = 1-((1-mean(CFDataVenusNorm))./0.44);

figure(3)
DosageDataEGFPNorm = 1-((1-CFDataEGFPNorm)./0.44);
DosageDataVenusNorm = 1-((1-CFDataVenusNorm)./0.44);
hold on
plot(1,DosageDataEGFPNorm,'go','MarkerFaceColor','g','MarkerEdgeColor','none','MarkerSize',10)
errorbar(1,mean(DosageDataEGFPNorm),std(DosageDataEGFPNorm),'ko','CapSize',0)
plot(2,DosageDataVenusNorm,'yo','MarkerFaceColor','y','MarkerEdgeColor','none','MarkerSize',10)
errorbar(2,mean(DosageDataVenusNorm),std(DosageDataVenusNorm),'ko','CapSize',0)
plot([0.5 2.5],[1 1],'k--')
%xlim([0.5 2.5])
xticks([1,2])
xticklabels({'EGFP-Bcd','Venus-Bcd'})
ylabel('protein dosage normalized to wt')

% calculate the dosage with respect to each other
VenusDosage./EGFPDosage


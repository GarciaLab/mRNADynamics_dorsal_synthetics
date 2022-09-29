function exploreBasicModel(c,nInactive,nOff,piEntry,tCycle,KDs)

% exploreBasicModel(1.5,0,5,0,8,logspace(3.2,5,8))
%fiveOffSteps(linspace(0,3800,18),1.5,KDs(1),nInactive,nOff,piEntry,[],tCycle,[])

% theta = [c,kd,nInactive,nOff,piEntry,piExit,tcycle]




piExit = 0;
dorsalVals = linspace(0,3800,18);
dorsalMidBin = dorsalVals(1:end) + diff(dorsalVals(1:2));
%KDs = logspace(3,5,7);
FractionPalette = cbrewer('seq', 'PuBu', length(KDs));
OnsetPalette = cbrewer('seq', 'OrRd', length(KDs));

% figFrac = figure;
% axFrac = axes(figFrac);
% figOn = figure;
% axOn = axes(figOn);
% 
% hold(axFrac,'on')
% hold(axOn,'on')
% figure
% hold on= [];
figure
Onsets = [];
for i = 1:length(KDs)
    
    kd = KDs(i);
    fractionColor = FractionPalette(i,:);
    onsetColor = OnsetPalette(i,:);
    
    
    theta = [c,kd,nInactive,nOff,piEntry,piExit,tCycle];
    modelOpts.TimeVariantDorsalValues = [];
    modelOpts.nSims = 1;
    fraction_onset = BasicModel_masterEq_DorsalTrace_AR(dorsalVals,theta, modelOpts);
    Fractions(i,:) = fraction_onset(:,1);
    Onsets(i,:) = fraction_onset(:,2);
    
    yyaxis left
    plot([0 dorsalMidBin],[zeros(i,1) Fractions],'-','Color',fractionColor,'LineWidth',2)
    ylim([0 1])
    ylabel('fraction')

    yyaxis right
    plot(dorsalMidBin,Onsets,'-','Color',onsetColor,'LineWidth',2)
    ylim([0 8])
    ylabel('onset')
    
end




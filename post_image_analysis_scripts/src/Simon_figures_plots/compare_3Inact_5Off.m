function compare_3Inact_5Off(c,piEntry,tCycle,KDs)


modelOpts.nSims=1;
modelOpts.TimeVariantDorsalValues=[];

dorsalVals = linspace(0,3800,18);
midBindorsalValues = dorsalVals(1:end-1) + diff(dorsalVals(1:2));
%c=1;
nInactive = 3;
nOff = 5;
%piEntry = 1;
%tCycle = 8;
%KDs = logspace(3.2,5,8);

theta = [c,KDs(1),nInactive,nOff,piEntry,0,tCycle];
% theta = [c,kd,nInactive,nOff,piEntry,piExit,tcycle]

fraction_onset = BasicModel_masterEq_DorsalTrace_AR(dorsalVals,theta, modelOpts);
ThreeInactStepsFiveOffSteps(dorsalVals,theta,modelOpts)

hold on
yyaxis left
plot(midBindorsalValues,fraction_onset(:,1),'bo')
yyaxis right
plot(midBindorsalValues,fraction_onset(:,2),'ro')
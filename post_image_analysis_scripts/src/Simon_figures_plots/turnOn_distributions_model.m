
c = 0.55;
kd = 1300;
TotalTime = 7.1;
NOffStates = 3;
theta = [c, kd, nan, NOffStates, nan,nan, TotalTime];
dorsalVals = linspace(0,3500,20);

fraction_onset = BasicModel_masterEq(dorsalVals,theta,[]);

function fraction_onset2 = ThreeInactStepsFiveOffSteps(dorsalVals,theta,options)
% theta = [c,kd,nInactive,nOff,piEntry,piExit,tcycle]


c = theta(1);
kd = theta(2);
nInactive = theta(3);
nOff = theta(4);
piEntry = theta(5);
piExit = theta(6);
tCycle = theta(7);


midBindorsalValues = dorsalVals(1:end-1) + diff(dorsalVals(1:2));
totalCycleTime = 12; %minutes
timeStepsPerMin = 10;
timeArraySize = totalCycleTime.*timeStepsPerMin;
t = linspace(0,totalCycleTime,timeArraySize);
k1single = (c*(dorsalVals./kd)./(1 + dorsalVals./kd));
endOfCycle = tCycle.*timeStepsPerMin;

%% calculate solution

k1 = repmat(k1single,numel(t), 1)'; % repeat the KD array so that it matches the dimensions of the time array
k2 = repmat(piEntry,numel(t), 1)';

EndState = 1+(1/24).*exp(1).^((-1).*(k1+k2).*t).*(k1+(-1).*k2).^(-7).*(12.* ... 
 exp(1).^(k2.*t).*k2.^5.*(2.*(21.*k1.^2+(-7).*k1.*k2+k2.^2)+2.*k1.* ...
  (6.*k1.^2+(-7).*k1.*k2+k2.^2).*t+k1.^2.*(k1+(-1).*k2).^2.*t.^2)+( ...
  -1).*exp(1).^(k1.*t).*k1.^3.*(24.*(k1.^4+(-7).*k1.^3.*k2+21.* ...  
k1.^2.*k2.^2+(-35).*k1.*k2.^3+35.*k2.^4)+24.*(k1+(-1).*k2).*k2.*( ...
  k1.^3+(-6).*k1.^2.*k2+15.*k1.*k2.^2+(-20).*k2.^3).*t+12.*(k1+(-1) ...
  .*k2).^2.*k2.^2.*(k1.^2+(-5).*k1.*k2+10.*k2.^2).*t.^2+4.*(k1+(-4) ...
  .*k2).*(k1+(-1).*k2).^3.*k2.^3.*t.^3+(k1+(-1).*k2).^4.*k2.^4.* ...
 t.^4));

% get the fraction of nuclei that are in the last state by the end of the
% available time
fraction_onset2(:,1) = EndState(:,endOfCycle);

% calculate the expected value of the time to reach the last (ON) state
Truncatedy6 = EndState(:,1:endOfCycle+1);
fraction_onset2(:,2) = sum(diff(Truncatedy6,1,2).*t(1:endOfCycle),2)./sum(diff(Truncatedy6,1,2),2);

% yyaxis left
% plot([0 midBindorsalValues],[0;fraction_onset2(:,1)],'b')
% ylim([0 1])
% 
% yyaxis right
% plot(midBindorsalValues,fraction_onset2(:,2),'r')
% ylim([0 8])

yyaxis left
plot(dorsalVals,fraction_onset2(:,1),'b')
ylim([0 1])

yyaxis right
plot(dorsalVals,fraction_onset2(:,2),'r')
ylim([0 8])



   
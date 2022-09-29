function fraction_onset2 = fiveOffSteps(dorsalVals,c,kd,n_entry,n_off,piEntry,piExit,tCycle,options)

midBindorsalValues = dorsalVals(1:end-1) + diff(dorsalVals(1:2));
totalCycleTime = 12; %minutes
timeStepsPerMin = 10;
timeArraySize = totalCycleTime.*timeStepsPerMin;
t = linspace(0,totalCycleTime,timeArraySize);
k = (c*(dorsalVals./kd)./(1 + dorsalVals./kd));
endOfCycle = tCycle.*timeStepsPerMin;

% for i = 1:length(k)
%     k1=k(i);
%     
% %     y1 = exp(1).^((-1).*k1.*t);
% %     y2 = exp(1).^((-1).*k1.*t).*k1.*t;
% %     y3 = (1/2).*exp(1).^((-1).*k1.*t).*k1.^2.*t.^2;
% %     y4 = (1/6).*exp(1).^((-1).*k1.*t).*k1.^3.*t.^3;
% %     y5 = (1/24).*exp(1).^((-1).*k1.*t).*k1.^4.*t.^4;
% 
%     y6 = 1+(1/24).*exp(1).^ ((-1).*k1.*t) .* ((-24)+(-1).*k1.*t.*(24+k1.*t.*(12+k1.*t.*(4+k1.*t))));
%     %y6 = 100*y6;
%     fraction_onset(i,2) = sum(diff(y6(1:endOfCycle+1)).*t(1:endOfCycle))/sum(diff(y6(1:endOfCycle+1)));
%     fraction_onset(i,1) = y6(endOfCycle);
%     
% %     tot = y1+y2+y3+y4+y5+(y6/100); % sanity check
% end

%% no looping
% now try getting rid of the loop over KDs and calculate everything at
% once.

k0 = repmat(k,numel(t), 1)'; % repeat the KD array so that it matches the dimensions of the time array
y6State = 1+(1/24) .* exp(-k0.*t) .* ((-24)-1*k0.*t.*(24+k0.*t.*(12+k0.*t.*(4+k0.*t))));

% get the fraction of nuclei that are in the last state by the end of the
% available time
fraction_onset2(:,1) = y6State(:,endOfCycle);

% calculate the expected value of the time to reach the last (ON) state
Truncatedy6 = y6State(:,1:endOfCycle+1);
fraction_onset2(:,2) = sum(diff(Truncatedy6,1,2).*t(1:endOfCycle),2)./sum(diff(Truncatedy6,1,2),2);

figure
hold on
yyaxis left
plot(dorsalVals,fraction_onset2(:,1),'b')
ylim([0 1])

yyaxis right
plot(dorsalVals,fraction_onset2(:,2),'r')
ylim([0 8])

% yyaxis left
% plot([0 midBindorsalValues],[0;fraction_onset2(:,1)],'b')
% ylim([0 1])
% 
% yyaxis right
% plot(midBindorsalValues,fraction_onset2(:,2),'r')
% ylim([0 8])

   
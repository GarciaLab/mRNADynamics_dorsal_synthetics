

t = 0:10:8*60;

% three parallel reversible steps, a, b and c.
% each step has their own forward (k1) and reverse (k2) rate constants

k1a = 0.01;
k2a = 0;% 0.001;
k1b = 0.005;
k2b = 0;%0.001;
k1c = 0.01;
k2c = 0;%0.001;


OnA_t = (k1a-(exp((-k1a-k2a).*t).*k1a))./(k1a+k2a);
OnB_t = (k1b-(exp((-k1b-k2b).*t).*k1b))./(k1b+k2b);
OnC_t = (k1c-(exp((-k1c-k2c).*t).*k1c))./(k1c+k2c);
All3On = OnA_t .* OnB_t .* OnC_t;

%plot(t,All3On)

DlRange = logspace(1,3.5,20);

figure(1)
hold on
DlKd = 100000;
counter = 1;
palette = viridis(length(DlRange));
for D = DlRange
    Color = palette(counter,:);
    K1 = DlBinding(D,DlKd,k1);
    k1a = K1;
    OnA_t = (k1a-(exp((-k1a-k2a).*t).*k1a))./(k1a+k2a);
    OnA_t = 1-exp(-k1a.*t);
    OnB_t = 1-exp(-k1b.*t);
    OnC_t = 1-exp(-k1c.*t);
    
    All3On = OnA_t .* OnB_t .* OnC_t;
    pOn(counter) = All3On(end);
    
    AvgTimeOn(counter) = sum(diff(All3On).*t(1:end-1))/sum(diff(All3On)); %expected value

    plot(t./60,All3On,'.-k')
    plot(t./60,OnA_t,'-b')
    plot(t./60,OnB_t,'-g')
    plot(t./60,OnC_t,'-y')

    %plot([AvgTimeOn(counter) AvgTimeOn(counter)]./60,[0 1],'Color',Color)
    
    counter = counter+1;
    
end
hold off
xlabel('time into NC12')
ylabel('fraction nuclei in the ON state')

%
figure(2)
plot(DlRange,pOn,'r')
xlabel('[Dl]')
ylabel('fraction of activer nuclei')
ylim([0 1.1])

figure(3)
plot(DlRange,AvgTimeOn./60,'b')
xlabel('[Dl]')
ylabel('mean turn on time (min)')
ylim([0 8.2])
    

    
    
    


function K1 = DlBinding(DlConcentration,DlKd,k1)
    Occupancy = (DlConcentration/DlKd)./(1+(DlConcentration/DlKd));
    K1 = k1*Occupancy;
end

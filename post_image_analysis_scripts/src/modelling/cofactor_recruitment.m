%%
Kd = 100; % dorsal binding affinity in AUs
D = logspace(-0.7,4,500) ; % Dorsal concentration in units of Kd

Kp = 40000; % RNAP binding affinity in AUs
P = 10 ; % RNAP concentration in units of Kp

Kc = 10000; % cofactor binding affinity 
C = 1 ; % Cofactor concentration in units of Kc

w1 = 20000;
w2 = 10000;

r = 1; % transcription rate when the promoter is active

Num = ((P/Kp) + ((C/Kc)*(P/Kp)*w2) + (D./Kd)*(C/Kc)*(P/Kp)*w1*w2)
Den = Num + 1 + (D./Kd) + (C/Kc) + (D./Kd)*(C/Kc)*w1
R = Num ./ Den;

figure(1)
plot(D,R)
set(gca,'XScale','log')

figure(2)
hold on
for Kd = logspace(2,5,7)
    Num = ((P/Kp) + ((C/Kc)*(P/Kp)*w2) + (D./Kd)*(C/Kc)*(P/Kp)*w1*w2)
    Den = Num + 1 + (D./Kd) + (C/Kc) + (D./Kd)*(C/Kc)*w1
    R = Num ./ Den;
plot(D,R)
end
hold off
set(gca,'XScale','log')


%%
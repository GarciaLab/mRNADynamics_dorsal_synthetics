function fraction_onset = fiveOffSteps(dorsalVals,c,kd,n_entry,n_off,piEntry,piExit,tCycle,options)

t = linspace(0,12,120);
k = (c*(dls./kd) ./ (1 + dls./kd));

for i = 1:length(k)
    k1=k(i);
    
%     y1 = exp(1).^((-1).*k1.*t);
%     y2 = exp(1).^((-1).*k1.*t).*k1.*t;
%     y3 = (1/2).*exp(1).^((-1).*k1.*t).*k1.^2.*t.^2;
%     y4 = (1/6).*exp(1).^((-1).*k1.*t).*k1.^3.*t.^3;
%     y5 = (1/24).*exp(1).^((-1).*k1.*t).*k1.^4.*t.^4;

    y6 = 1+(1/24).*exp(1).^ ((-1).*k1.*t) .* ((-24)+(-1).*k1.*t.*(24+k1.*t.*(12+k1.*t.*(4+k1.*t))));
    y6 = 100*y6;
    fraction_onset(i,1) = sum(diff(y6(1:81)).*t(1:80))/sum(diff(y6(1:81)));
    fraction_onset(i,2) = y6(80)./100;
    
%     tot = y1+y2+y3+y4+y5+(y6/100)
end

   
   
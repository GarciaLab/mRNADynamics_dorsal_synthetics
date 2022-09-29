function accumulatedmRNA = tfDrivenAccumulatedmRNA(c, kd, R, d)

%c- rate increase per unit dorsal
%d- dorsal concentration
%t- time into nc12
%kd- kd of dorsal in simple binding mode
%R- strength of transcriptional response given Dorsal is bound
t = 10; %this t represents a cutoff time of 10 minutes into nc12. 

 accumulatedmRNA = (1/24).*d.*(d+kd).^(-1).*R.*(24.*t+exp(-c.*d.*(d+kd).^( ...
  -1).*t).*(120.*c.^(-1).*d.^(-1).*(d+kd)+96.*t+36.*c.*d.*(d+kd).^( ...
  -1).*t.^2+8.*c.^2.*d.^2.*(d+kd).^(-2).*t.^3+c.^3.*d.^3.*(d+kd).^( ...
  -3).*t.^4));
 
end
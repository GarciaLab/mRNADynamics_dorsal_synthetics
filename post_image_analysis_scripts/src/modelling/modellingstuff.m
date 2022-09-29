 %this function models the curve of affinity vs. accumulated mrna
%inputs: duration of cycle. (mean production rate dmrna/dt at each time
%point and dv bin). total number of nuclei within each dv bin. total number
%of nuclei transcribing within each dv bin. 

deltaT_s = 10; %s
frame = 0:1:60; 
time_s = frame * deltaT_s;

frame_on = 18;
dorsal = 0:250:4500;  %aus 
KDs = 5:1:5000;
fluo_offset = 0;
N_nuclei = 18; %in nc12

%constant hill mRNA rate
hill = dorsalFitFunction('hill');
n = 1; %hill coefficient
scale = 1;

for i = 1:length(KDs)

    parameters_hill = [scale, KDs(i), n, fluo_offset];
    mRNA_rate_const = hill(parameters_hill, dorsal);
%     plot(dorsal, mRNA_rate_const)

    mRNA_rate = zeros(length(frame), length(dorsal)); 

    %toy scenario where transcription starts ~3 minutes into the cycle, lasts
    %til the end and is constant.
    for k = 1:length(dorsal)
        mRNA_rate(frame_on:end, k) = mRNA_rate_const(k);
    end


%     [X,Y] = meshgrid(dorsal, time_s);
%     surf(X,Y,mRNA_rate) 
%     ylabel('time (s)');
%     xlabel('[Dorsal] (au)');
%     zlabel('mRNA production rate')

    mRNA_rate_nozeros = mRNA_rate;
    mRNA_rate_nozeros(mRNA_rate_nozeros==0) = nan;

    mRNA_rate_mean = nanmean(mRNA_rate_nozeros, 1); %time averaged production rate
    mRNA_rate_mean( isnan(mRNA_rate_mean) ) = 0;
    mRNA_rate_mean_DVIntegral(i) = trapz(dorsal, mRNA_rate_mean); %DV integrated time averaged production rate
    
end

figure(1)
plot(KDs, mRNA_rate_mean_DVIntegral);
set(gca, 'XDir','reverse')
set(gca, 'XScale', 'log');
xlabel('log KD (aus)')
ylabel('mean accumulated fluorescence');

%%


%constant hill mRNA rate
hill = dorsalFitFunction('hill');
n = 1; %hill coefficient
scale = 1;
KD = 500;

parameters_hill = [scale, KD, n, fluo_offset];
mRNA_rate_const = hill(parameters_hill, dorsal);

figure(2)
plot(dorsal, mRNA_rate_const)
xlabel('dorsal')
ylabel('mRNA production rate, hill')

mRNA_rate = zeros(length(frame), length(dorsal)); 

%toy scenario where transcription starts ~3 minutes into the cycle, lasts
%til the end and is constant.
for k = 1:length(dorsal)
    mRNA_rate(frame_on:end, k) = mRNA_rate_const(k);
end

figure(3);
[X,Y] = meshgrid(dorsal, time_s);
surf(X,Y,mRNA_rate) 
ylabel('time (s)');
xlabel('[Dorsal] (au)');
zlabel('mRNA production rate')


%%

%4 state kinetic model

dlon = 1;
dlon_prime = 1;
polon = 1;
polon_prime = 1;
dloff = 1;
dloff_prime = 1;
poloff = 1;
poloff_prime = 1;
r_poldl = 1;
r_pol = 1;

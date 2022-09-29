
% close all force;


dmax = 5000;
nPlots = 10;
dl = 500;
dls = logspace(1, log10(dmax)); %aus
t = linspace(0, 10)'; %mins
kd = 500;
kds = logspace(2, 4, nPlots);
cs = logspace(-1, 2, nPlots);
pi1s = logspace(-2, 1, nPlots);

R = 500;
% c = .002;
c = 40;

t_cycle = 10; %min

nSteps = 7;
nSims = 1E3;
nOffStates = 5;
onstate = nOffStates+1;
silentstate = onstate+1;
occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );
pi0 = 1; %min-1
pi1 = .1; %min-1

dls = 700; %au

cmap = colormap(viridis(nPlots));

cs = 2;

noExit = false;

if noExit
    pi1s = 0;
    dls = logspace(1, log10(dmax)); %aus
    kds = logspace(2, 4, nPlots);
    cs = logspace(-1, 2, nPlots);
else
    dls = logspace(1, log10(dmax)); %aus
    kds = logspace(2, 4, nPlots);
    cs = logspace(-1, 2, nPlots);
    pi1s = logspace(-2, 1, nPlots);
end

mfpts = nan(length(dls), length(kds), length(pi1s));
factive = nan(length(dls), length(kds), length(pi1s));

for m = 1:length(cs)
    for k = 1:length(pi1s)
        for i = 1:length(dls)
            for j = 1:length(kds)
                
                pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
                pi1 = pi1s(k); %min-1
                
                [fpt_on_observed, factive_temp] = averagePaths(nSims, nSteps, pi0, pi1, onstate, silentstate, t_cycle);
                
                mfpts(i, j, k, m) = mean(fpt_on_observed);
                factive(i,j,k,m) = factive_temp;
                
            end
        end
    end
end

dt = mfpts(:, 10, :, :) - mfpts(:, 4, :, :); %kds(10)=10k, kds(4)=400

dt = repmat(dt, [1 length(kds) 1 1]);

[~, dropboxfolder] = getDorsalFolders;

if noExit
    save([dropboxfolder, 'tfdriven_paramsearch.mat']);
    load([dropboxfolder, 'tfdriven_paramsearch.mat'])
else
    save([dropboxfolder, 'tfdrivenexit_paramsearch.mat'])
    load([dropboxfolder, 'tfdrivenexit_paramsearch.mat'])   
end

try
    plotTFDrivenParams(factive, dt, mfpts, 'nPoints', 2E4)
catch
    plotTFDrivenParams(factive, dt, mfpts);
end

figure;
tiledlayout('flow')
for j = 1:length(cs)
    nexttile;
    for k = 1:length(pi1s)
        %     plot(kds, mfpts(40, :, 1));
        plot(kds, mfpts(1, :, k, j), 'LineWidth', 2, 'Color', cmap(k, :));
        hold on
    end
    xlabel('K_D (au)')
    ylabel('mean time to turn on (min)')
    %     leg =legend(num2str(round2(pi1s')));
    %     title(leg, '\pi_1');
    %     set(gca, 'XScale', 'log')
    title(['c = ', num2str(cs(j))])
end

figure;
plot(pi1s, mfpts(:, 5));
xlabel('super off transition rate (min-1)')
ylabel('mean time to turn on (min)')

%
% figure;
% plot(dls, mfpts(:, 5));
% xlabel('[Dl] (au)')
% ylabel('mean time to turn on (min)')

function [states, times] = makePath(nSteps, pi0, pi1, onstate, silentstate,  tau_on, tau_exit, ind, states, times)

state = 1;

for step = 2:nSteps
    
    if ind == 1 && state < onstate
        
        state =  state + 1;
        tau = tau_on(step);
        states(step) = state;
        times(step) = times(step-1) + tau;
        
    elseif ind == 2 || state == onstate
        
        states(step) = silentstate;
        tau = tau_exit(step);
        times(step) = times(step-1) + tau;
        return;
        
    end
    
end

end

function [fpt_on_observed, factive] = averagePaths(nSims, nSteps, pi0, pi1, onstate, silentstate, t_cycle)

fpt_on = nan(1, nSims);

%     fpt_on_observed = [];
%     duration = [];

taus = [exprnd(pi0^-1, [1, nSteps, nSims]) %on
    exprnd(pi1^-1, [1, nSteps, nSims]) %exit
    ];

tau_on = squeeze(taus(1, :, :));
tau_exit = squeeze(taus(2, :, :));

[~, ind] = min([taus(1,:, :); taus(2,:, :)]);


states = [1, zeros(1, nSteps-1)];
times = zeros(1, nSteps);

for k = 1:nSims
    
    [states, times] = makePath(nSteps, pi0, pi1, onstate, silentstate, tau_on(:, k), tau_exit(:, k), ind(:, k), states, times);
    
    ton = times(states==onstate);
    if ~isempty(ton)
        fpt_on(k) = ton;
    end
    
end

fpt_on_observed = fpt_on(fpt_on < t_cycle);
factive = length(fpt_on_observed) / nSims;

%     duration = t_cycle-fpt_on_observed;

%     histogram(fpt_on_observed, 'Normalization', 'pdf');
%     legend(['<\tau>=', num2str(round2(mean(fpt_on_observed))), ' min']);
%     xlabel('\tau (min)')
%     ylabel('probability density function')


end


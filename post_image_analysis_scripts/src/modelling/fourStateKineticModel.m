function fourStateKineticModel()

close all;

dl = 10:10:500;

figure;
tiledlayout;

nexttile;
dmrnadt = [];
for k = 1:length(dl)
    [t, y] = solveFour(dl(k));
    dmrnadt(k,:) = y(:, 5);
    plot(t, dmrnadt(k, :));
    hold on
end
xlabel('time (min')
ylabel('dmRNA/dt (au/min)')

nexttile;
plot(dl, dmrnadt(:, 60));
xlabel('[Dl] (nM)')
ylabel('dmRNA/dt at t = 6 mins (au/min)')





end

function [t, y] = solveFour(dl)

    %pol and dorsal are assumed to have the same binding & unbinding
    %behaviors and concentrations. kons are calculated with berg purcell 
    %4pi*.3*60 formula multiplied to dorsal concentration
    if nargin < 1
     dl = 200; %nm
    end
    
    pol = 20; %nm?
    
    %we'll bound everything between 1 and 10^5 min^-1; 
    k01= 230*pol; %pol kon
    k02= 230*dl; %dl on
    k10= 20; %pol off. assuming same as dl
    k13= 230*pol; %dl on with pol on
    k20= 20; %dl off. should be right order of magnitude
    k23= 230*pol; %pol on with dl on
    k31= 20; %dl off with pol on. assuming this = k20.
    k32= 20; %pol off with dl on
    
    %non-dimensionalize on rates by dividing by the average dorsal
    %concentration. then all of the above are in min^-1
    
    K = [(-k10-k20) k01 k02 0
            k10 (-k01-k31) 0 k13
            k20 0 (-k02-k32) k23
            0 k31 k32 (-k13-k23)];
    
    R = 500; %aus/min
    tspan = [0:.1:10];
    y0 = [1 0 0 0 0];
    
    [t,y] = ode45(@(t,y) odefcn(y(:), K, R), tspan, y0(:));
    
    displayFigures = false;
    if displayFigures
        
        figure;
        tiledlayout(1, 2);
        nexttile;
        plot(t, y(:, 5));
        set(gca, 'XScale', 'log');
        xlim([10^-2, 10]); %min
        ylabel('dmRNA/dt (au/min)')
        xlabel('time (min)')
        
        nexttile;
        plot(t,y(:,1), 'r',t,y(:,2), 'k', t, y(:, 3), 'b', t, y(:, 4), 'g')
            legend('state 1', 'state 2', 'state 3', 'state 4');
        ylabel('prob. of being in state')
        xlabel('time (min')
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
    end
    
end

function dydt = odefcn(y, K, R)

    dydt(1:4, 1) = K*y(1:4);
    dydt(5, 1) = R * ( y(2)+y(3) ./ ( y(1) + y(2) + y(3) + y(4) ) );
    
end
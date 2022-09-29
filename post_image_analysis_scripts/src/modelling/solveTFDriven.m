function [t, y, dmrnadt, mrna] = solveTFDriven(K, R, displayFigures)
    
    %this can be made into an anonymous function to pass to fitting
    %functions (lsqnonlin, etc.) like this-
    % fun = @(K, R) solveFour(K, R, false)
    
    %pol and dorsal are assumed to have the same binding & unbinding
    %behaviors and concentrations. kons are calculated with berg purcell 
    %4pi*.3*60 formula multiplied to dorsal concentration
    
    tspan = 0:.001:10;
    y0 = [1 0 0 0];
    
    opts = odeset('Vectorized','on', 'Jacobian', K);

    [t,y] = ode23s(@(t,y) K*y, tspan, y0(:), opts);
    
    %invoking the occupancy hypothesis
    dmrnadt = R * ( y(:, 1)+y(:, 4) ./ ( y(:, 1) + y(:, 2) + y(:, 3) + y(:, 4) ) );
    mrna = cumtrapz(t, dmrnadt);
    
    
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
            legend('state 0', 'state 1', 'state 2', 'state 3');
        ylabel('prob. of being in state')
        xlabel('time (min')
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
    end
    
end
function fourStateKineticModel()

% close all;

R = 200; %aus/min
    
dl_au = 10:50:3000; %au
co_au = 500; %averageish [Dl]
dl = dl_au/co_au; %dimensionless

dmrnadts = [];
mrnas = [];
t0 = [];
for k = 1:length(dl)
    
    pol = 20;
    %we'll bound everything between 1 and 10^5 min^-1; 
    k01= 3.77*pol; %pol kon
    k23= 3.77*pol; %pol on with dl on
    k02= 3.77*dl(k); %dl on
    k13= 3.77*dl(k); %dl on with pol on
    k10= 20; %pol off
    k32= .1; %pol off with dl on
    k20= 20; %dl off
    k31= 20; %dl off with pol on

    %non-dimensionalize on rates by dividing by the average dorsal
    %concentration. then all of the above are in min^-1

    K = [(-k10-k20) k01 k02 0
            k10 (-k01-k31) 0 k13
            k20 0 (-k02-k32) k23
            0 k31 k32 (-k13-k23)];
    
    [t, y, dmrnadt, mrna] = solveFour(K, R, false);
    
    t0(k, :) = t; %#ok<*AGROW>
    mrnas(k, :) = mrna;
    dmrnadts(k, :) = dmrnadt;
    state0(k, :) = y(:, 1);
    state1(k, :) = y(:, 2);
    state2(k, :) = y(:, 3);
    state3(k, :) = y(:, 4);
    
end

figure;
tiledlayout('flow');

nexttile;
palette = parula(length(dl));
for k = 1:length(dl)
    plot(t0(k, :), dmrnadts(k, :), 'Color', palette(k,:));
    hold on
end
xlabel('time (min)')
ylabel('dmRNA/dt (au/min)')
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
cbh = colorbar;
set(cbh,'YTickLabel',dl_au)
ylabel(cbh, '[Dl] (au)');

nexttile;
palette = parula(size(dmrnadts, 2));
t00 = t0(1, :);
for k = 1:size(dmrnadts, 2)
    plot(dl_au, dmrnadts(:, k), 'Color', palette(k,:));
    hold on
end
xlabel('[Dl] (au)')
ylabel('dmRNA/dt (au/min)')
c = colorbar;
set(c,'YTickLabel',t00)
ylabel(c, 'time (min)');



nexttile;
palette = parula(size(mrnas, 2));
t00 = t0(1, :);
for k = 1:size(mrnas, 2)
    plot(dl_au, mrnas(:, k),'Color', palette(k,:));
    hold on
end
xlabel('[Dl] (au)')
ylabel('accumulated mRNA (au)')
c = colorbar;
set(c,'YTickLabel',t00)
ylabel(c, 'time (min)');



nexttile;
plot(dl_au, state0(:, 60), 'r', dl_au, state1(:, 60), 'k', dl_au, state2(:, 60), 'b', dl_au, state3(:, 60), 'g')
xlabel('[Dl] (au)')
ylabel('prob. of state at t = 6 mins (au/min)')
legend('state 0', 'state 1', 'state 2', 'state 3');
set(gca, 'YScale', 'log');

% b = uicontrol('Parent',f,'Style','slider',...
%               'value',20, 'min',0, 'max',100, 'Callback', @update_k1);


end
% 
% function update_k1(source, ~)
%     
% end
function plotSimpleWeak()

close all;

dmax = 5000;
d = [0:.1:dmax];
f = @(x, p) p{1}.*(((x./p{2}).*p{3})./(1+ (x./p{2})+ ((x./p{2}).*p{3})));

R = 500;
KD = 500;
w = 3;

nPlots = 7;
cmap = colormap(viridis(nPlots));

figure; tiledlayout('flow')

nexttile;
ws = logspace(0, 1, nPlots);
for k = 1:nPlots
    p = {R, KD, ws(k)};
    plot(d, f(d, p), 'LineWidth', 2, 'Color', cmap(k, :));
    hold on
end
set(gca, 'XScale', 'log')
leg = legend(gca, num2str(round(ws(:), 2, 'significant')), 'Location','northwest' );
title(leg, '\omega''')
xlabel('[Dorsal] (au)')
ylabel('R (au/min)')
xlim([0, dmax]);
ylim([0, R*1.1]);



nexttile;
KDs = logspace(2, 4, nPlots);
for k = 1:nPlots
    p = {R, KDs(k), w};
    plot(d, f(d, p), 'LineWidth', 2, 'Color', cmap(k, :));
    hold on
end
set(gca, 'XScale', 'log')
leg = legend(gca, num2str(round(KDs(:), 2, 'significant')), 'Location','northwest' );
title(leg, 'K_D')
xlabel('[Dorsal] (au)')
ylabel('R (au/min)')
xlim([0, dmax]);
ylim([0, R*1.1]);



nexttile;
Rs = logspace(1, 3, nPlots);
for k = 1:nPlots
    p = {Rs(k), KD, w};
    plot(d, f(d, p), 'LineWidth', 2, 'Color', cmap(k, :));
    hold on
end
set(gca, 'XScale', 'log')
leg = legend(gca, num2str(round(Rs(:), 2, 'significant')), 'Location','northwest' );
title(leg, 'R_{max}')
xlabel('[Dorsal] (au)')
ylabel('R (au/min)')
xlim([0, dmax]);
ylim([0, R*1.1]);



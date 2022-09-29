
kds = logspace(1, 6);
gs = [];
vwkds = [];
dgkds = [];
for k = 1:length(kds)
    disp(kds(k))
    [results,chain,s2chain] = fitstuff_mcmc2glob('nSimu', 1E5, 'metric',...
        'fluo', 'minKD', 0, 'maxKD', kds(k), 'displayFigures', false, 'wb', false);
    y = chainstats(chain, results);
    gs(k) = mean(y(:, end));
    vwkds(k) = y(8, 1);
    dgkds(k) = y(2, 1);
end

figure; tiledlayout('flow');
nexttile;
plot(kds, gs)
xlabel('kd upper limit (au)')
ylabel('mean geweke')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

nexttile;
plot(nsims, vwkds)
ylabel('predicted 1DgVW KD (au)')
xlabel('kd upper limit (au)')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
nexttile;
plot(nsims, dgkds)
ylabel('predicted 1Dg KD (au)')
xlabel('kd upper limit (au)')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

figure; tiledlayout('flow');
nsims = round(logspace(2, 6));
gs2 = [];
vwkds2 = [];
dgkds2 = [];
for k = 1:length(nsims)
    [results,chain,s2chain] = fitstuff_mcmc2glob('nSimu', nsims(k), 'metric',...
        'fluo', 'minKD', 0, 'maxKD', 1E4, 'displayFigures', false, 'wb', false);
    y = chainstats(chain, results);
    gs2(k) = mean(y(:, end));
    vwkds2(k) = y(8, 1);
    dgkds2(k) = y(2, 1);
end

nexttile;
plot(nsims, gs2)
xlabel('nsims')
ylabel('mean geweke')
nexttile;
plot(nsims, vwkds2)
ylabel('predicted 1DgVW KD (au)')
xlabel('nsims')
nexttile;
plot(nsims, dgkds2)
ylabel('predicted 1Dg KD (au)')
xlabel('nsims')
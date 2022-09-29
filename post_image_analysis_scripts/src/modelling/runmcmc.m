close all force;

nSims = 1E5; 
minw = 1E-2;
maxw = 1E1;

% [results, chain, s2chain] =fitstuff_mcmc2glob('nSimu', nSims, 'expmnt', 'affinities', 'metric', 'fluo', 'minKD', 0,...
%     'maxKD', 9000, 'noOff', true, 'minw', minw, 'maxw', maxw, 'minR', 1E-3, 'maxR', 1E4, 'displayFigures', false);
% 
% kd_1dg = results.mean(2)
% 
% [results, chain, s2chain] =fitstuff_mcmc2glob('nSimu', 1E5, 'expmnt', 'phases', 'metric', 'fraction',...
%     'noOff', true, 'minw', minw, 'maxw', maxw, 'minR', 1E-3, 'maxR', 1E4, 'fixedKD', kd_1dg);
% 
% 
% [results, chain, s2chain] =fitstuff_mcmc2glob('nSimu', 1E5, 'expmnt', 'phases', 'metric', 'fluo',...
%     'noOff', true, 'minw', minw, 'maxw', maxw, 'minR', 1E-3, 'maxR', 1E4, 'fixedKD', kd_1dg);

%%

[results, chain, s2chain] =fitstuff_mcmc2glob('nSimu', nSims, 'expmnt', 'affinities', 'metric', 'fluo', 'minKD', 0,...
    'maxKD', 9000, 'noOff', true, 'minw', minw, 'maxw', maxw, 'minR', 1E-3, 'maxR', 1E4, 'displayFigures', true,...
    'enhancerSubset', {'1Dg11', '1DgS2', '1DgW'}, 'scoreSubset', [6.23, 5.81, 5.39]);

fixedR = results.mean(end);
fixedw = results.mean(1);

[results, chain, s2chain] =fitstuff_mcmc2glob('nSimu', nSims, 'expmnt', 'affinities', 'metric', 'fluo', 'minKD', 0,...
    'maxKD', 9000,  'noOff', true, 'minw', minw, 'maxw', maxw, 'minR', 1E-3, 'maxR', 1E4, 'displayFigures', true,...
    'enhancerSubset', {'1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'}, 'scoreSubset', [5.13, 4.80, 4.73, 4.29], 'fixedR',fixedR, 'fixedw', fixedw );

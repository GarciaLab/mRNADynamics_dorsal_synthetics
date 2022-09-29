function fitstuff_mcmc2()

close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};



%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = nan(length(enhancers), 2);
xrange(1,:)  = [1000, 2250];  %1Dg
xrange(2, :) = [500, 1750]; %upper limit for %1DgS
xrange(3, 1) = 500; %lower limit for 1DgW
xrange(4, 2) = 1500; %upper limit for %1DgAW


nSets = length(enhancers);
xo = {};
yo = {};
xs = {};
ys = {};
dsid = [];
T = [];
Y = [];
for k = 1:nSets
    xo{k} = dorsalResultsDatabase.dorsalFluoBins( ...
        strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}) );
    yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo( ...
        strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    T = [T; xs{k}];
    Y = [Y; ys{k}];
end

data.ydata = [xs{1}, ys{1} ];


%%
% Refine the first guess for the parameters with fminseacrh and calculate residual variance as an estimate of the model error variance.
modd = standardizeModelForGramm(dorsalFitFunction('hill'));
model2 = dorsalFitFunction('hill');
x = xs{1};
y = ys{1};
%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
k00 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 1; -y_max];
ub = [y_max*2; Inf; 8; y_max*10];



[mdl,gof,~] =fit(x,y,modd,...
    'StartPoint',k00, 'Lower',lb, 'Upper',ub);
k0 = coeffvalues(mdl)';
mse = (gof.rmse).^2;
%
params = cell(1, 3);
for i = 1:length(k0)
    params{1, i} = {['k', num2str(i)],k0(i), 0};
end

%set max n of hill function to 10
params{1, 3}{4} = 10;

model = struct;
model.ssfun = @funss;
model.sigma2 = mse;

options.nsimu = 10000;
options.updatesigma = 1;
[results,chain,s2chain] = mcmcrun(model,data,params,options);

figure(1); clf
mcmcplot(chain,[],results,'chainpanel')
subplot(2,2,4)
mcmcplot(sqrt(s2chain),[],[],'dens',2)
title('error std')

% Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.

chainstats(chain,results)

figure(2); clf
xx = (0:1:1E4)';
yy = model2(mean(chain), xx);
plot(xo{1},yo{1},'s',xx,yy,'-')
xlim([min(xx), max(xx)])
ylim([min(yy), max(yy)]);
legend({'obs','fit'},'Location','best')
title('Data and fitted model')


%ideally, these guys look like ellipses. if certain parameters give weird
%shapes, it might mean those parameters should be removed from the model if
%possible
figure(3); clf
mcmcplot(chain,[],results,'pairs', .5);

% commented out because these graphs are inside the pairs figure
% figure(4); clf
% mcmcplot(chain,[],results,'denspanel',2);

end

function ss = funss(k, data)
% sum-of-squares
time = data.ydata(:,1);
Aobs = data.ydata(:,2);
m = dorsalFitFunction('hill');
Amodel = m(k, time);

ss = sum((Aobs-Amodel).^2);
end
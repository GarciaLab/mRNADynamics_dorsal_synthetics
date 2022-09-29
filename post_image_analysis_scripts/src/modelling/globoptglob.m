function globoptglob(varargin)


[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')

expmnt = "affinities"; %affinities, phases or phaff
md = "simpleweak"; %simpleweak, simpleweakdimer, repression, tfdriven, artifact, fourstate
metric = "fluo"; %fraction, fluo
lsq = false;
noOff = true;
nSimu = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
minKD = 200;
maxKD = 1E4;
minw = 1E-2; %1E-2
maxw = 1E1; %1E2
minR = 10;
maxR = 1E3;
displayFigures = true;
wb = true;
fixedKD = NaN; %if this value isn't nan, this value will determine the fixed KD parameter and KD won't be fitted
fixedOffset = NaN;
fixedR = NaN;
fixedw = NaN;
enhancerSubset = {};
scoreSubset = [];
positionSubset = [];
useBatches = false; %fit all the data across embryos and not just means 

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if expmnt == "phaff"
    lsq = false;
    noOff = true;
end

enhancers_1dg = {'1Dg11'};
enhancers_aff =  {'1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
enhancers_ph = {'1Dg-5', '1Dg-8D'};
scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
positions = [0, -5, -8]';

if strcmpi(expmnt, 'affinities')
    enhancers = [enhancers_1dg, enhancers_aff];
elseif strcmpi(expmnt, 'phases')
    enhancers = [enhancers_1dg, enhancers_ph];
elseif expmnt=="phaff"
    enhancers = [enhancers_1dg, enhancers_aff, enhancers_ph];
end

if ~isempty(enhancerSubset)
    enhancers = enhancerSubset;
    scores= scoreSubset;
    positions = positionSubset;
end


%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = getXRange(enhancers, expmnt);

nSets = length(enhancers);
xo = {};
yo = {};
xs = {};
ys = {};
dsid = [];
T = [];
Y = [];
T_batch = [];
Y_batch = [];
dsid_batch = [];
for k = 1:nSets
    cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k});
    xo{k} = dorsalResultsDatabase.dorsalFluoBins(cond);
    if metric == "fraction"
        yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
        yo_batch{k} = dorsalResultsDatabase.fracFluoEmbryo(cond, :);
    elseif metric == "fluo"
        yo{k} = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond);
        yo_batch{k} = dorsalResultsDatabase.allMaxFluoEmbryo(cond, :);
    end
    
     xo_batch{k} = repmat(xo{k}, [1, size(yo_batch{k}, 2)]);
    [xs_batch{k}, ys_batch{k}]= processVecs(xo_batch{k}, yo_batch{k}, xrange(k, :));
    
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    dsid_batch = [dsid_batch; k*ones(size(xo{k}))];

    if isnan(xrange(k, 1))
        xrange(k, 1) = min(xo{k});
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = max(xo{k});
    end
    
    assert(~isempty(xs{k}));
    
    T = [T; xs{k}];
    Y = [Y; ys{k}];
    T_batch = [T_batch; xo{k}];
    Y_batch = [Y_batch; ys_batch{k}];
end

data.ydata = [T, Y];
data.dsid = dsid;
data.X =  [T dsid];

data_batch = {};
for k = 1:size(Y_batch, 2)
    data_batch{k}.ydata = [T_batch Y_batch(:, k)];
    data_batch{k}.X =  [T_batch dsid_batch];
    data_batch{k}.dsid = dsid_batch;
end

if useBatches
    disp('Doing batched fits');
    data = data_batch;
end

%%


%rate, kd, hill, y offset
y_max = nanmax(Y(:));
x_max = max(T);

[p0, lb, ub, names] = getInits(expmnt, md, metric ,x_max, y_max, nSets,...
    'minKD', minKD, 'maxKD', maxKD, 'minw', minw, 'maxw', maxw, 'minR', minR, 'maxR', maxR, 'subset', ~isempty(enhancerSubset));


if lsq
    % Refine the first guess for the parameters with fminseacrh and calculate residual variance as an estimate of the model error variance.
    [k0, mse] = globfit2('expmnt',expmnt, 'metric', metric, 'md', md, 'maxKD', maxKD, 'displayFigures', displayFigures);
else
    k0 = p0;
end


if noOff
    %ignore the offset in the initial parameters and bounds
    p0(names=="offset") = [];
    lb(names=="offset") = [];
    ub(names=="offset") = [];
    k0(names=="offset") = [];
    names(names=="offset") = [];
end

simpleWeakOptions = struct('noOff', noOff, 'fraction', metric=="fraction",...
    'dimer', contains(md, "dimer"), 'expmnt', expmnt);
mdl = @(x, p) simpleweak(x, p, simpleWeakOptions);

rng(1,'twister'); %set the rng seed so we get the same results every run of this function

%%
x = [];
%fminunc
ssfun = @(params) sum( (data.ydata(:,2)-mdl(data.X(:,1), params)).^2 );
rf2 = ssfun;
x0 = p0;
% [x(:, 1),ff,flf,of] = fminunc(rf2,x0);


%patternsearch
[x(:, 1),fp,flp,op] = patternsearch(rf2,x0, [], [], [], [], lb, ub);

%genetic algorithm
nvars = length(p0);
% initpop = 10*randn(20,nvars) + repmat(x0,[1,20])';
% opts = optimoptions('ga','InitialPopulationMatrix',initpop);
% [x(:, 3),fga,flga,oga] = ga(rf2,nvars,[],[],[],[],[],[],[],opts);
[x(:, 2),fga,flga,oga] = ga(rf2,nvars,[],[],[],[], lb, ub);


%particle swarm
% opts = optimoptions('particleswarm','InitialSwarmMatrix',initpop);
% [x(:, 4),fpso,flgpso,opso] = particleswarm(rf2,nvars,[],[],opts);
[x(:, 3),fpso,flgpso,opso] = particleswarm(rf2,nvars,lb, ub);

%surrogate optimization
opts = optimoptions('surrogateopt','PlotFcn',[]);
[x(:, 4),fsur,flgsur,osur] = surrogateopt(rf2,lb,ub,opts);


%global search
problem = createOptimProblem('fmincon','objective',rf2,...
    'x0',p0, 'lb',lb,'ub',ub);
gs = GlobalSearch;
[x(:, 5),fg,flg,og] = run(gs,problem);

figure; tiledlayout('flow')
%fminunc is pretty lame, so let's exclude it from the histogram. 
% x = x(:, 2:end);
for k = 1:nvars
    nexttile;
    if contains(names(k), 'KD')
        bw = 100;
    elseif contains(names(k), 'w')
        bw = .1;
    elseif contains(names(k), 'R')
        bw = 20;
    end
    histogram(x(k, :), 'BinWidth', bw);
    title(names(k));
end


figure;
nAlgos = 5;
for k = 1:nvars
    if k ~=1
     yyaxis left
    else
        yyaxis right
    end
    plot(k, x(k, :), '.', 'MarkerSize', 30);
    hold on
%     title(names(k));
end
xlim([0, nvars+1])
% set(gca, 'YScale', 'log')

1


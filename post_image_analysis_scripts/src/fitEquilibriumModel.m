function [results, chain, s2chain, data, modelOpts] = fitEquilibriumModel(varargin)

wb = true;
nSteps = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
fixKD = false;
batchedAffinities = true;
preRun = false;
bin = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

[~, dropboxfolder] = getDorsalFolders;

%1Dg data (highest affinity)
datPath = dropboxfolder + "\manuscript\window\basic\dataForFitting\archive\";
load(datPath + "binMidValues.mat", "binMidValues");
load(datPath + "MaxFluosPerEmbryoAll.mat", "MaxFluosPerEmbryoAll");


if batchedAffinities
    
    
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
    data_batch = {};
    
    
    for k = 1:numel(enhancers)
        
        
        %needs to be nobs x ny
        MaxFluosPerEmbryo = MaxFluosPerEmbryoAll{k};
        
        if bin
            X = binMidValues;
            F = nanmean(MaxFluosPerEmbryo, 1);
        else
            X = repmat(binMidValues, 1, size(MaxFluosPerEmbryo, 1));
            F = MaxFluosPerEmbryo(:)';
            X(isnan(F)) = [];
            F(isnan(F)) = [];
            X(F==0) = [];
            F(F==0) = [];
        end
        data{k}.ydata = [X; F]';
    end
    
else
    MaxFluosPerEmbryo = MaxFluosPerEmbryoAll{1};
    %needs to be nobs x ny
    X = repmat(binMidValues, 1, max(size(MaxFluosPerEmbryo)));
    F = MaxFluosPerEmbryo(:)';
    X(isnan(F)) = [];
    F(isnan(F)) = [];
    data.ydata = [X; F]';
end




%%
rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
% options.drscale = 2; % a high value (5) is important for multimodal parameter spaces.
options.drscale = 5;
options.waitbar = wb; %the waitbar is rate limiting sometimes
options.nsimu = nSteps; %should be between 1E3 and 1E6
options.updatesigma = 1; %honestly don't know what this does
% options.method = 'mh';
%

names = ["R", "kd", "w"];
p0 = [100, 1E3, 1];
lb = [1, 1E2, 1E-2];
ub = [1E4, 1E5, 1E2];

fixedKDs = [208, 155, 977, 872, 872, 1487, 3694];

params = cell(1, length(p0));
pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance

for k = 1:length(names)
    
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    localflag = 0; %is this local to this dataset or shared amongst batches?
    
    
    if fixKD && names(k) == "kd"
        targetflag = 0;
    end
    
    
    %we allow each enhancer to have a different kd but share every other
    %param
    if batchedAffinities && names(k) == "kd"
        localflag = 1;
    end
    
    if fixKD && batchedAffinities && names(k) == "kd"
        params{1, k}= {names(k), fixedKDs, lb(k), ub(k), pri_mu, pri_sig, targetflag, localflag};
    else
        params{1, k}= {names(k), p0(k), lb(k), ub(k), pri_mu, pri_sig, targetflag, localflag};
    end
    
end


modelOpts = struct;


model = struct;

%mdl inputs are 1. dorsal 2. {amp, KD, omegaDP, offset}. Any of the elements
%of 2 are allowed to be column vectors.
mdl =  @(x, p) p(1).*(((x./p(2)).*p(3))./(1+ (x./p(2))+ ((x./p(2)).*p(3))));

model.modelfun   = @(x, p) p(1).*(((x.ydata(:, 1)./p(2)).*p(3))./(1+ (x.ydata(:, 1)./p(2))+ ((x.ydata(:, 1)./p(2)).*p(3)))); %use mcmcrun generated ssfun

results = [];
% [results,~,~,~]=mcmcrun(model,data,params,options,results);
if preRun
    [results,~,~,~]=mcmcrun(model,data,params,options,results);
end
[results,chain,s2chain,~]=mcmcrun(model,data,params,options,results);

burnInTime = .25; %let's burn the first 25% of the chain just in case
chainLen = size(chain, 1);
chain = chain(round(burnInTime*chainLen):chainLen, :);
if ~isempty(s2chain)
    s2chain = s2chain(round(.25*chainLen):chainLen, :);
end



%%
% close all force;

chainfig = figure(); clf
mcmcplot(chain,[],results,'chainpanel')

% Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.
%geweke is a measure of whether things converged between 0 and 1.
chainstats(chain,results)


try
    %ideally, these guys look like ellipses. if certain parameters give weird
    %shapes, it might mean those parameters should be removed from the model if
    %possible
    pairFig = figure; clf
    % mcmcplot(chain,[],results,'pairs', .5);
    mcmcplot(chain,[],results, 'pairs', 4);
    
    %Make the scatter plot more visible
    ax = gca;
    ax.Children(4).MarkerSize = 4;
    ax.Children(4).Color = [.4 .4 1];
end
%
% figure;
out = mcmcpred(results,chain,[],data, results.modelfun);


%%
lims = [0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995];
figure;
tiledlayout(1, numel(enhancers));

for k = 1:length(out.predlims)
    nn = (size(out.predlims{k}{1},1) + 1) / 2;
    plimi = out.predlims{k};
    x{k} = data{k}.ydata(:, 1);
    
    nexttile;
    for j = 1:length(plimi)
     
        %         y_lower{k, j} = plimi{j}(4,:); %for 25-75%
%         y_upper{k, j}  = plimi{j}(6,:); %for 25-75 %
        y_lower{k, j} = plimi{j}(1,:); %for 95%
        y_upper{k, j}  = plimi{j}(2*nn-1,:); %for 95%
        
        y_median{k, j} = plimi{j}(nn,:);
        
        yl{k, j} = nan(1, length(binMidValues));
        yu{k, j} = nan(1, length(binMidValues));
        ym{k, j} = nan(1, length(binMidValues));
        
        for i = 1:length(binMidValues)
            yl{k, j}(i) = nanmedian(y_lower{k, j}(x{k}==binMidValues(i)));
            yu{k, j}(i) = nanmedian(y_upper{k, j}(x{k}==binMidValues(i)));
            ym{k, j}(i) = nanmedian(y_median{k, j}(x{k}==binMidValues(i)));
            
            temp_f = MaxFluosPerEmbryoAll{k}(:, i);
            temp_f(temp_f == 0) = nan;
            nFVals{k}(i) = sqrt(length(temp_f(~isnan(temp_f))));
            
        end
        
        yl2 =  yl{k, j};
        yu2 = yu{k, j};
        x7 = binMidValues(~isnan(yl2));
        yl2(isnan(yl2)) = [];
        yu2(isnan(yu2)) = [];
        fillyy(x7, yl2, yu2,[0.9 0.9 0.9]);
        hold on
        nFVals{k}(nFVals{k}==0) = nan;
        yyy = MaxFluosPerEmbryoAll{k};
        yyy(yyy==0) = nan;
        y4 = nanmean(yyy, 1);
        ye4 =  nanstd(yyy, 1)./nFVals{k};
        x4 = binMidValues(~isnan(y4));
        y4(isnan(y4)) = [];
        ye4(isnan(ye4)) = [];
        errorbar(x4, y4, ye4,...
            'o', 'CapSize', 0)
        ylabel('max. fluo.(a.u.)')
        xlabel('Dorsal concentration (a.u.)')
        ylim([0, 450]);
        
        y5 = ym{k, j};
        x5 = binMidValues(~isnan(y5));
        y5 = y5(~isnan(y5));
        plot(x5, y5, '-')
        
        axis square;
    end
    
end
%%




% mcmcpredplot(out);

%%
if isempty(results.mean)
    results.mean = mean(chain);
    warning('Results mean was empty. Using chain mean.')
end

theta_mean = getAllThetas(names, results);


if fixKD && batchedAffinities
    for k = 1:length(enhancers)
        yy{k} = mdl(binMidValues, [theta_mean(1), fixedKDs(k), theta_mean(3)]);
    end
else
    if ~iscell(theta_mean)
        yy = mdl(binMidValues, theta_mean);
    else
        for k = 1:length(theta_mean)
            yy{k} = mdl(binMidValues, theta_mean{k});
        end
    end
end


if ~iscell(yy)
    figure;
    tiledlayout('flow')
    nexttile;
    
    nFVals = sqrt(length(MaxFluosPerEmbryo(~isnan(MaxFluosPerEmbryo))));
    errorbar(binMidValues, nanmean(MaxFluosPerEmbryo, 1), nanstd(MaxFluosPerEmbryo, 1)/nFVals)
    hold on
    plot(binMidValues, yy(:, 1))
    
    
    % nexttile
    % scatter(sims.factive(:), sims.mfpts(:))
    % scatter(yy(:, 1), yy(:, 2))
    % xlim([0, 1.1])
    
    % nexttile;
    % plot(data.ydata(1, :), data.ydata(2, :), data.ydata(2, :))
    % ylabel('fraction')
    % xlabel('dl')
    % zlabel('onset')
else
    figure;
    tiledlayout('flow')
    for k = 1:length(yy)
        nexttile;
        yyy = MaxFluosPerEmbryoAll{k};
        yyy(yyy==0) = nan;
        nFVals = sqrt(length(yyy(~isnan(yyy))));
        errorbar(binMidValues, nanmean(yyy, 1), nanstd(yyy./nFVals, 1))
        hold on
        plot(binMidValues, yy{k}, 'LineWidth', 2)
        
        ylabel('factive')
        xlabel('dl')
        %         legend('data', 'sim')
        title(enhancers{k})
        
        
    end
end



save(dropboxfolder + "\simulations\" +"equilibrium"+"_"+datestr8601+".mat")

disp('debug stop')
keyboard;
end
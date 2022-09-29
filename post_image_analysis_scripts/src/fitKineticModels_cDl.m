function [results, chain, s2chain, data, modelOpts] = fitKineticModels_cDl(varargin)

% To do:


wb = true;
nSteps = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
nSims = 1E3; %number of simulations for the kinetic barrier model (not the number of mcmc walker steps).
exitOnlyDuringOffStates = true; %determines connectivity of the markov graph
modelType = "entryexit"; %choices- entryexit, entry, exit, basic
fun= "table"; %also 'sim', 'imhomo', 'master', 'masterInhomo'
variableStateNumber = false;
fixKD = false;
batchedAffinities = false;
fixTCycle = false;
preRun = false;
bin = false;
piForm = "cdl";

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
load(datPath + "FractionsPerEmbryoAll.mat", "FractionsPerEmbryoAll");
load(datPath + "OnsetsPerEmbryoAll.mat", "OnsetsPerEmbryoAll");


if batchedAffinities
    
    %this setting being true wouldn't make sense, so let's fix it.
    fixKD = false;
    
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
    data_batch = {};
    
    
    for k = 1:numel(enhancers)
        
        
        %needs to be nobs x ny
        FractionsPerEmbryo = FractionsPerEmbryoAll{k};
        OnsetsPerEmbryo = OnsetsPerEmbryoAll{k};
        
        if bin
            X = binMidValues;
            F = nanmean(FractionsPerEmbryo, 1);
            T =  nanmean(OnsetsPerEmbryo, 1);
        else
            X = repmat(binMidValues, 1, max(size(FractionsPerEmbryo)));
            F = FractionsPerEmbryo(:)';
            T = OnsetsPerEmbryo(:)';
            X(isnan(F)) = [];
            T(isnan(F)) = [];
            F(isnan(F)) = [];
        end
        data{k}.ydata = [X; F; T]';
    end
    
else
    
    FractionsPerEmbryo = FractionsPerEmbryoAll{1};
    OnsetsPerEmbryo = OnsetsPerEmbryoAll{1};
    %needs to be nobs x ny
    X = repmat(binMidValues, 1, max(size(FractionsPerEmbryo)));
    F = FractionsPerEmbryo(:)';
    T = OnsetsPerEmbryo(:)';
    X(isnan(F)) = [];
    T(isnan(F)) = [];
    F(isnan(F)) = [];
    data.ydata = [X; F; T]';
end



%%
saveStr = modelType;
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
if variableStateNumber
    saveStr = saveStr + "variableStateNumber";
end
% sims = load(dropboxfolder +  "\simulations\archive\" + "tf_paramsearch_"+saveStr+"_.mat", 'params', 'factive', 'mfpts');
if fun == "table"
    sims = load(dropboxfolder +  "\simulations\" + "tf_paramsearch_"+saveStr+"_.mat", 'params', 'factive', 'mfpts');
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
%names = ["c", "kd" , "nentrystates", "moffstates", "pentry", "pexit", "tcycle"];
names = ["c", "nentrystates", "moffstates", "pentry", "pexit", "tcycle"];
if modelType == "entryexit"
    if ~batchedAffinities
        p0 = [10, 5, 5, 1, 1, 6.8];
        lb = [1E-1,  0, 1, 1E-1, 1E-2, 5];
        ub = [1E2,  12, 12, 1E3, 1E1, 10];
    else
        p0 = [1,  1, 1, 1, .01, 8];
        lb = [1E-1, 0, 1, 1E0, 1E-3, 7];
        ub = [1E2,  12, 12, 1E2, 1E-2, 9];
    end
elseif modelType == "entry"
    p0 = [10,  5, 5, 1, 0, 8];
    lb = [1E-2,0, 1, 1E-1, 0, 5];%pentry lower than .1 causes crash
    ub = [1E2,  12, 12, 1E1, 0, 10];
elseif modelType == "basic"
    p0 = [.001, 0, 4, 1E10, 0, 7.1];
    lb = [1E-5,  0, 1, 1E10, 0, 4];
    ub = [1E1, 0, 12, 1E10, 0, 9];
end

if variableStateNumber
    lb(1) = 10;
end

params = cell(1, length(p0));
pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance

for k = 1:length(names)
    
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    localflag = 0; %is this local to this dataset or shared amongst batches?
    
    if  contains(names(k), "states") && ~variableStateNumber
        targetflag = 0;
    end
    
    if ~contains(modelType, "exit")
        if names(k) == "pexit"
            targetflag = 0;
            p0(k) = 0;
        end
    end
    
    if ~contains(modelType, "entry")
        if names(k) == "pentry" || names(k) == "nentrystates"
            targetflag = 0;
            p0(k) = 0;
        end
    end
        
    if fixTCycle && names(k) == "tcycle"
        targetflag = 0;
    end
    
    %we allow each enhancer to have a different c but share every other
    %param
    if batchedAffinities && names(k) == "c"
        localflag = 1;
    end
    
    params{1, k}= {names(k), p0(k), lb(k), ub(k), pri_mu, pri_sig, targetflag, localflag};
end


modelOpts = struct;
modelOpts.modelType = modelType;
modelOpts.nSims = nSims;

if fun == "table"
    modelOpts.sims = sims;
elseif ~contains(fun,"master")
    modelOpts.exitOnlyDuringOffStates = true;
    
    gpurng(1, "ThreeFry"); %fastest gpu rng
    
    nSilentStates = contains(modelType, 'exit');
    if variableStateNumber
        nentries = ub(3);
        moffs = ub(4);
    else
        nentries = p0(3);
        moffs = p0(4);
    end
    nStates = nentries + moffs + 1 + nSilentStates;
    n_dls = length(binMidValues);
    modelOpts.r_vec = gpuArray(rand(1, nentries*nSims +...
        ( (moffs+1) * nSims * n_dls ) +...
        nSilentStates*((nStates-1) * nSims), 'single'));
end

model = struct;

if fun == "table"
    mdl = @(x, p) kineticFunForFits_table(x, p, modelOpts);
elseif fun== "sim"
    mdl = @(x, p) kineticFunForFits_sim_vec_gpu_customrnd(x, p, modelOpts);
elseif fun == "inhomo"
    mdl = @(x, p)  timesim_interp_alldl(x, p, modelOpts);
elseif fun == "master" && modelType == "basic"
    mdl = @(x, p)  BasicModel_masterEq(x, p, modelOpts);
elseif fun=="masterInhomo" && modelType == "basic"
    fullMatFileName = datPath + '/DorsalFluoTraces.mat';
    load(fullMatFileName);
    modelOpts.TimeVariantDorsalValues = [DorsalFluoTraces.meanDorsalFluo];
    modelOpts.TimeVariantAbsoluteTimes = DorsalFluoTraces(1).absoluteTime; %in seconds
    modelOpts.middleBinValues = [DorsalFluoTraces.binValue];
    if piForm == "cdl"
        modelOpts.piForm = "cdl";
    end
    mdl = @(x, p)  BasicModel_masterEq_DorsalTrace_AR(x, p, modelOpts);
end
% mdl = @(x, p) kineticFunForFits_sim(x, p, modelOpts);
model.modelfun   = mdl;  %use mcmcrun generated ssfun

if ~batchedAffinities
    model.ssfun = @(theta, data) sum( (data.ydata(:, 2:end) -...
        mdl(data.ydata(:, 1), theta) ).^2, 'omitnan');
end

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
out = mcmcpred(results,chain,[],data, mdl);


%%
lims = [0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995];
figure;
for k = 1:length(out.predlims)
    nn = (size(out.predlims{k}{1},1) + 1) / 2;
    plimi = out.predlims{k};
    x{k} = data{k}.ydata(:, 1);
    
    nexttile;
    for j = 1:length(plimi)
        
        %         y_lower{k, j} = plimi{j}(1,:);
        y_lower{k, j} = plimi{j}(4,:);
        %  y_upper{k, j}  = plimi{j}(2*nn-1,:);
        y_upper{k, j}  = plimi{j}(6,:);
        
        y_median{k, j} = plimi{j}(nn,:);
        
        yl{k, j} = nan(1, length(binMidValues));
        yu{k, j} = nan(1, length(binMidValues));
        ym{k, j} = nan(1, length(binMidValues));
        
        for i = 1:length(binMidValues)
            yl{k, j}(i) = nanmedian(y_lower{k, j}(x{k}==binMidValues(i)));
            yu{k, j}(i) = nanmedian(y_upper{k, j}(x{k}==binMidValues(i)));
            ym{k, j}(i) = nanmedian(y_median{k, j}(x{k}==binMidValues(i)));
            
            temp_f = FractionsPerEmbryoAll{k}(:, i);
            nFVals{k}(i) = sqrt(length(temp_f(~isnan(temp_f))));
            temp_o = OnsetsPerEmbryoAll{k}(:, i);
            nOVals{k}(i) = sqrt(length(temp_o(~isnan(temp_o))));
        end
        if j == 1
            yyaxis left
        elseif j ==2
            yyaxis right
        end
        hold on
        fillyy(binMidValues,yl{k, j},yu{k, j},[0.9 0.9 0.9]);
        hold on
        plot(binMidValues, ym{k, j})
        hold on
        if j ==1
            hold on
%             nFVals = sqrt(length(FractionsPerEmbryoAll{k}(~isnan(FractionsPerEmbryoAll{k}))));
            nFVals{k}(nFVals{k}==0) = nan;

            errorbar(binMidValues, nanmean(FractionsPerEmbryoAll{k}, 1), nanstd(FractionsPerEmbryoAll{k}./nFVals{k}, 1),...
                'o', 'CapSize', 0)
            ylabel('factive')
            xlabel('dl')
            ylim([0, 1.1])
        elseif j == 2
            hold on
%             nOVals = sqrt(length(OnsetsPerEmbryoAll{k}(~isnan(OnsetsPerEmbryoAll{k}))));
            nOVals{k}(nOVals{k}==0) = nan;

            errorbar(binMidValues, nanmean(OnsetsPerEmbryoAll{k}, 1), nanstd(OnsetsPerEmbryoAll{k}, 1)./nOVals{k},...
                'o', 'CapSize', 0)
            ylabel('onset')
            xlabel('dl')
            ylim([0, 8.5])
            title(enhancers{k})

        end
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


if ~iscell(theta_mean)
    yy = results.modelfun(binMidValues, theta_mean);
else
    for k = 1:length(theta_mean)
        yy{k} = results.modelfun(binMidValues, theta_mean{k});
    end
end

if ~iscell(yy)
    figure;
    tiledlayout('flow')
    nexttile;
    nFVals = sqrt(length(FractionsPerEmbryo(~isnan(FractionsPerEmbryo))));
    errorbar(binMidValues, nanmean(FractionsPerEmbryo, 1), nanstd(FractionsPerEmbryo, 1)/nFVals)
    hold on
    plot(binMidValues, yy(:, 1))
    
    ylabel('factive')
    xlabel('dl')
    legend('data', 'sim')
    nexttile;
    nOVals = sqrt(length(OnsetsPerEmbryo(~isnan(OnsetsPerEmbryo))));
    errorbar(binMidValues, nanmean(OnsetsPerEmbryo, 1), nanstd(OnsetsPerEmbryo, 1)/nOVals)
    
    hold on
    plot(binMidValues, yy(:, 2))
    
    ylabel('onset')
    xlabel('dl')
    legend('data', 'sim')
    
    % nexttile
    % scatter(sims.factive(:), sims.mfpts(:))
    % scatter(yy(:, 1), yy(:, 2))
    % xlim([0, 1.1])
    % ylim([0, 8.1])
    
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
        yyaxis left
        nFVals = sqrt(length(FractionsPerEmbryoAll{k}(~isnan(FractionsPerEmbryoAll{k}))));
        errorbar(binMidValues, nanmean(FractionsPerEmbryoAll{k}, 1), nanstd(FractionsPerEmbryoAll{k}./nFVals, 1))
        hold on
        plot(binMidValues, yy{k}(:, 1))
        ylim([0, 1.1])
        
        ylabel('factive')
        xlabel('dl')
        %         legend('data', 'sim')
        title(enhancers{k})
        yyaxis right
        nOVals = sqrt(length(OnsetsPerEmbryoAll{k}(~isnan(OnsetsPerEmbryoAll{k}))));
        errorbar(binMidValues, nanmean(OnsetsPerEmbryoAll{k}, 1), nanstd(OnsetsPerEmbryoAll{k}, 1)./nOVals)
        
        hold on
        plot(binMidValues, yy{k}(:, 2))
        
        ylabel('onset')
        xlabel('dl')
        %         legend('data', 'sim')
        title(enhancers{k})
        
    end
end



% %%
% try
%     n = size(chain, 1);
%     kd = results.theta(2)*ones(n, 1);
%     nentries = results.theta(3)*ones(n, 1);
%     moffs = results.theta(4)*ones(n, 1);
%     piexits = results.theta(6)*ones(n, 1);
%     theta = horzcat(chain(:, 1), kd, nentries, moffs, chain(:, 2), piexits);
%
%     yyy = nan(length(binMidValues),length(results.mean),n);
%     for k = 1:n
%         yyy(:, :, k) = results.modelfun(binMidValues, theta(k, :));
%     end
% end
%
% nexttile;
% yyy_f = squeeze(yyy(:, 1, :));
% yyy_o = squeeze(yyy(:, 2, :));
%
% scatter(yyy_f(:), yyy_o(:));
% xlim([0, 1.1])
% ylim([0, 8.2])
%
% nexttile;
% % scatter(sims.factive(:), sims.mfpts(:))
% xlim([0, 1.1])
% ylim([0, 8.2])

% %%
% figure;
% a = [];
% for k = 1:length(OnsetsPerEmbryoAll)
%     a = [a; OnsetsPerEmbryoAll{k}(:)];
%     hold on
% end
% histogram(a(:));
%

save(dropboxfolder + "\simulations\" +modelType+fun+nSteps+"_"+datestr8601+".mat")

disp('debug stop')
keyboard;
end
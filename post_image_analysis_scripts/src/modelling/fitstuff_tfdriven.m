%% Fraction active
close all;

% enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3'};

scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73]';

    [~, resultsFolder] = getDorsalFolders;
    load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')


%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = nan(length(enhancers), 2);
xrange(1,:)  = [1000, 2250];  %1Dg
xrange(2, :) = [500, 1750]; %upper limit for %1DgS
xrange(3, 1) = 500; %lower limit for 1DgW
xrange(4, 2) = 1500; %upper limit for %1DgAW

%%
figure;
t = tiledlayout(1, length(enhancers), 'TileSpacing','Compact');
fits = struct;

for k = 1:length(enhancers)

    x = dorsalResultsDatabase.dorsalFluoBins( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}) );
    y = dorsalResultsDatabase.meanFracFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );
    y_error = dorsalResultsDatabase.seFracFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );
    
       model = @(c, KD, x) tfDrivenFractionActive(c, KD, x);

   %c kd 
   p0 = [.5, 500];
    lb = [0, 0]; 
   ub = [Inf, Inf];
  
    
    %remove nans from the data for fitting
    y_error(isnan(y)) = [];
    x(isnan(y)) = [];
    y(isnan(y)) = [];

    %interpolate the data to be able to fit and get CIs. otherwise the fit is
    %impossible since #params ~ #data points. 
    scale_interp = 2;
    xq = (min(x): mean(diff(x))/scale_interp : max(x) )';
    yq = interp1(x,y,xq);

    %restrict the range of the fit to a monotonically increasing region
    
    if isnan(xrange(k, 1))
        xrange(k, 1) = xq(1);
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = xq(end);
    end
    
    x1_ind = find(xq==xrange(k, 1));
    x2_ind = find(xq==xrange(k, 2));
    
    xin = xq(x1_ind:x2_ind);
    yin =  yq(x1_ind:x2_ind);
    mdl=fit(xin,yin,model,...
                'StartPoint',p0, 'Lower',lb, 'Upper',ub);

  % plotting
  
    ax = nexttile;

    errorbar(x, y, y_error);
    hold on
    plot(mdl)

    ci95 = predint(mdl,xin,0.95,'functional','off');
    plot(xin,ci95,'m--')

    ylim([0, 1]);
    xlim([0, 3250]);
    box(gca, 'on')
    set(gca, 'XLabel', []);
    set(gca, 'YLabel', []);
    
    if k~=length(enhancers)
        legend('hide')
    else
        legend({'Data \pm SEM','Fitted curve', 'Prediction intervals'},...
           'FontSize',8,'Location','northeast')
    end
    
    [coeffText, fits] = getModelText(mdl, fits, k);
    
    
    titleText = [ enhancers{k}; ['fit range: ', num2str(xq(x1_ind)), ' to ', num2str(xq(x2_ind)) ]; coeffText];
    title(ax, titleText);
    
end
%%

%assuming model form is the same for every graph so grabbing the last
%formula used
f = formula(mdl);
f = replace(f, [".", "(x)", "KD", "x"], ["", "x", "K_D", "[Dl]"]);
t.Title.String = "model: " + f;
xlabel(t,'Dorsal concentration (au)')
ylabel(t,'fraction active')


t = plotModelCoeffs(fits, scores, enhancers);
title(t, 'fraction active hill parameters')

%%

%% Accumulated fluorescence

figure;
t = tiledlayout(1, length(enhancers), 'TileSpacing','Compact');

for k = 1:length(enhancers)

    x = dorsalResultsDatabase.dorsalFluoBins( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}) );
    y = dorsalResultsDatabase.meanallmrnasEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );
    y_error = dorsalResultsDatabase.seallmrnasEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );

    model = @(c, KD, R, x) tfDrivenAccumulatedmRNA(c, KD, R, x);

    %remove nans from the data for fitting
    y_error(isnan(y)) = [];
    x(isnan(y)) = [];
    y(isnan(y)) = [];

    %interpolate the data to be able to fit and get CIs. otherwise the fit is
    %impossible since #params ~ #data points. 
    scale_interp = 2;
    xq = (min(x): mean(diff(x))/scale_interp : max(x) )';
    yq = interp1(x,y,xq);

    %restrict the range of the fit to a monotonically increasing region
    
    if isnan(xrange(k, 1))
        xrange(k, 1) = xq(1);
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = xq(end);
    end
    
    x1_ind = find(xq==xrange(k, 1));
    x2_ind = find(xq==xrange(k, 2));
    
    xin = xq(x1_ind:x2_ind);
    yin =  yq(x1_ind:x2_ind);
    
    
       %c kd R 
   p0 = [.5, 500, 500];
   lb = [0.01, 10, 1]; 
   ub = [2, 1E4, 1E3];

    
    mdl=fit(xin,yin,model,...
                'StartPoint',p0, 'Lower',lb, 'Upper',ub);

  % plotting
  
    ax = nexttile;

    errorbar(x, y, y_error);
    hold on
    plot(mdl)

    ci95 = predint(mdl,xin,0.95,'functional','off');
    plot(xin,ci95,'m--')

    ylim([0, 1250]);
    xlim([0, 3250]);
    box(gca, 'on')
    set(gca, 'XLabel', []);
    set(gca, 'YLabel', []);
    
    if k~=length(enhancers)
        legend('hide')
    else
        legend({'Data \pm SEM','Fitted curve', 'Prediction intervals'},...
           'FontSize',8,'Location','northeast')
    end
    
   [coeffText, fits] = getModelText(mdl, fits, k);
    
    titleText = [ enhancers{k}; ['fit range: ', num2str(xq(x1_ind)), ' to ', num2str(xq(x2_ind)) ]; coeffText];
    title(ax, titleText);
    
end

%assuming model form is the same for every graph so grabbing the last
%formula used
f = formula(mdl);
f = replace(f, [".", "(x)", "KD", "x"], ["", "x", "K_D", "[Dl]"]);
t.Title.String = "model: " + f;
xlabel(t,'Dorsal concentration (au)')
ylabel(t,'accumulated fluorescence (au)')


t = plotModelCoeffs(fits, scores, enhancers);
title(t, 'accumulated fluo hill parameters')

%% Fraction active
close all;

enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};

scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';

%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = nan(length(enhancers), 2);
xrange(1,:)  = [1000, 2250];  %1Dg
xrange(2, :) = [500, 1750]; %upper limit for %1DgS
xrange(3, 1) = 500; %lower limit for 1DgW
xrange(4, 2) = 1500; %upper limit for %1DgAW

model = standardizeModelForGramm(dorsalFitFunction('simpleWithPol'));


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
    
    %rate, kd, pkd, offset, omegadp
    y_max = nanmax(y(:));
    x_max = max(x);
    %specific for fraction active fitting
    p0 = [y_max; x_max/2 ; 1 ;0; 1];
    lb = [0; 1000; 0; 0; 1];
    ub = [min(1, y_max*2); Inf; Inf; y_max*10; Inf];
    
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


%% Max Fluorescence

figure;
t = tiledlayout(1, length(enhancers), 'TileSpacing','Compact');

for k = 1:length(enhancers)

    x = dorsalResultsDatabase.dorsalFluoBins( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}) );
    y = dorsalResultsDatabase.meanAllMaxFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );
    y_error = dorsalResultsDatabase.seAllMaxFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );


   %rate, kd, pkd, offset, omegadp
    y_max = nanmax(y(:));
    x_max = max(x);
    p0 = [y_max; x_max/2 ; 1; 1 ; 1];
    lb = [0; 1000; 0; 0; 1];
    ub = [y_max*2; Inf; Inf; y_max*10; Inf];

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

    ylim([0, 800]);
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
ylabel(t,'max fluorescence (au)')


%coefficient graphs
figure;
t = tiledlayout(1, length(fits(1).vals), 'TileSpacing','Compact');

for j = 1:length(fits(1).vals)
    nexttile;
    y = [];
    y_error_lower = [];
    y_error_upper = [];
    for k = 1:length(enhancers)
        y = [y, fits(k).vals(j)]; %#ok<AGROW>
        y_error_upper = [ y_error_upper, fits(k).cis(2, j) - fits(k).vals(j) ];
        y_error_lower = [ y_error_lower, fits(k).vals(j) - fits(k).cis(1, j) ];
    end
    plot(scores, y,'LineWidth', 2);
    hold on
    errorbar(scores, y, y_error_lower, y_error_upper);
    title(fits(1).names{j})
    ylim([.5*min(y), 1.5*max(y) ])
end

t = plotModelCoeffs(fits, scores, enhancers);
title(t, 'max fluo hill parameters')


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


    %rate, kd, pol/kd_pol, y offset, omega pol dl
    y_max = nanmax(y(:));
    x_max = max(x);
    p0 = [y_max; x_max/2 ; 1 ; 0; 1];
    lb = [0; 1000; 0; 0; 1];
    ub = [y_max*2; Inf; Inf; y_max*10; Inf];

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


%%
%let's just do this analysis for 1dg
f = standardizeModelForGramm(dorsalFitFunction('simpleWithPol'));

x = dorsalResultsDatabase.dorsalFluoBins( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, '1Dg11') );

scale_interp = 100;
x = (min(x): mean(diff(x))/scale_interp : max(x) )';      

R = 1;
KD = 2200;
PKDP = .0013;
offset = .013;
omegaDP = 1400;

KDs = 100:100:3100;

figure;
t = tiledlayout('flow');

nexttile;
palette = parula(length(KDs));
for k = 1:length(KDs)
    y = f(R, KDs(k), PKDP, offset, omegaDP, x);
    plot(x, y,  'Color', palette(k, :))
    hold on
end
c = colorbar;
set(c,'TickLabels', KDs(1:3:end))
ylabel(c, 'KD (au)');
hold off;
ylabel('normalized dmRNA/dt')
xlabel('[Dorsal] (au)');

omegaDPs = logspace(-1, 4);
nexttile;
palette = parula(length(omegaDPs));
for k = 1:length(omegaDPs)
    y = f(R, KD, PKDP, offset, omegaDPs(k), x);
    plot(x, y,  'Color', palette(k, :))
    hold on
end
c = colorbar;
set(c,'TickLabels', round(omegaDPs(1:5:end), 2, 'significant') )
ylabel(c, 'OmegaDP');
ylabel('normalized dmRNA/dt')
xlabel('[Dorsal] (au)');

PKDPs = logspace(-4, -2);
nexttile;
palette = parula(length(PKDPs));
for k = 1:length(PKDPs)
    y = f(R, KD, PKDPs(k), offset, omegaDP, x);
    plot(x, y,  'Color', palette(k, :))
    hold on
end
c = colorbar;
set(c,'TickLabels', round(PKDPs(1:5:end), 2, 'significant') )
ylabel(c, 'Pol/KDP');
ylabel('normalized dmRNA/dt')
xlabel('[Dorsal] (au)');


    %%
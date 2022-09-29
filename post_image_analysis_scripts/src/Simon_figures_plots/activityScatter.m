function activityScatter


%%
%this is where the data is
load('S:\Simon\Dropbox\DorsalSyntheticsDropbox\dorsalResultsDatabase.mat')

enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3'};
colorList = {'b','g','y','r','c','m'};

%% Fraction On (on X) vs max fluorescence (on Y)
figure
hold on
XforFit = [];
YforFit = [];
for k = 1:length(enhancers)
%     X = dorsalResultsDatabase.dorsalFluoBins( ...
%     strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
    % fraction on the X axis
    X = dorsalResultsDatabase.meanFracFluoEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    XError = dorsalResultsDatabase.seFracFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
       % max fluo Y axis
    Y = dorsalResultsDatabase.meanAllMaxFluoEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    YError = dorsalResultsDatabase.seAllMaxFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
    XforFit = [XforFit;X];
    YforFit = [YforFit;Y];
    
    errorbar(X,Y,YError,'o','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
    errorbar(X,Y,XError,'o','horizontal','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
   
end

% remove nans
X = XforFit(~isnan(XforFit));
Y = YforFit(~isnan(XforFit));

%fit to linear model
[p,S] = polyfit(X,Y,1);
f = polyval(p,linspace(0,1,20)); 
%plot(linspace(0,1,20),f,'k-')
% get R^2
[R,P,RL,RU] = corrcoef(X,Y);
Rsqrd = R(1,2);

%fit to line going through the origin
% slope = X\Y; % linear regression throug origin
% yfit = linspace(0,1,20)*slope; % calculate fitted line
mdl = fitlm(X,Y,'Intercept',false);
slope = mdl.Coefficients.Estimate;
yfit = linspace(0,1,20).*slope;
Rsqrd = mdl.Rsquared.Ordinary;
plot(linspace(0,1,20),yfit,'k-')

hold off
ylim([0 600])
xlim([0 1])
xlabel('fraction active')
ylabel('maximum fluorescence (AU)')
title(['R^2 = ' num2str(Rsqrd)])
%set(gca,'YScale','log'); set(gca,'XScale','log')





%% Fraction On (on X) vs accumulated fluo (on Y)
figure
hold on
XforFit = [];
YforFit = [];
for k = 1:length(enhancers)

    % fraction on the X axis
    X = dorsalResultsDatabase.meanFracFluoEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    XError = dorsalResultsDatabase.seFracFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
       % max fluo Y axis
    Y = dorsalResultsDatabase.meanallmrnasEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    YError = dorsalResultsDatabase.seallmrnasEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
    XforFit = [XforFit;X];
    YforFit = [YforFit;Y];
    
    errorbar(X,Y,YError,'o','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
    errorbar(X,Y,XError,'o','horizontal','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
   
end

% remove nans
X = XforFit(~isnan(XforFit));
Y = YforFit(~isnan(XforFit));

%fit to linear model
[p,S] = polyfit(X,Y,1);
f = polyval(p,linspace(0,1,20)); 
%plot(linspace(0,1,20),f,'k-')
% get R^2
[R,P,RL,RU] = corrcoef(X,Y);
Rsqrd = R(1,2);

%fit to line going through the origin
% slope = X\Y; % linear regression throug origin
% yfit = linspace(0,1,20)*slope; % calculate fitted line
mdl = fitlm(X,Y,'Intercept',false);
slope = mdl.Coefficients.Estimate;
yfit = linspace(0,1,20).*slope;
Rsqrd = mdl.Rsquared.Ordinary;
plot(linspace(0,1,20),yfit,'k-')

hold off
ylim([0 750])
xlim([0 1])
xlabel('fraction active')
ylabel('accumulated fluorescence (AU)')
title(['R^2 = ' num2str(Rsqrd)])
%set(gca,'YScale','log'); set(gca,'XScale','log')




%% max fluo (on X) vs accumulated fluo (on Y)
figure
hold on
XforFit = [];
YforFit = [];
for k = 1:length(enhancers)

    % fraction on the X axis
    X = dorsalResultsDatabase.meanAllMaxFluoEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    XError = dorsalResultsDatabase.seAllMaxFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
       % max fluo Y axis
    Y = dorsalResultsDatabase.meanallmrnasEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    YError = dorsalResultsDatabase.seallmrnasEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
    XforFit = [XforFit;X];
    YforFit = [YforFit;Y];
    
    errorbar(X,Y,YError,'o','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
    errorbar(X,Y,XError,'o','horizontal','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
   
end

% remove nans
X = XforFit(~isnan(XforFit));
Y = YforFit(~isnan(XforFit));

%fit to linear model
[p,S] = polyfit(X,Y,1);
f = polyval(p,linspace(0,1,20)); 
%plot(linspace(0,1,20),f,'k-')
% get R^2
[R,P,RL,RU] = corrcoef(X,Y);
Rsqrd = R(1,2);

%fit to line going through the origin
% slope = X\Y; % linear regression throug origin
% yfit = linspace(0,1,20)*slope; % calculate fitted line
mdl = fitlm(X,Y,'Intercept',false);
slope = mdl.Coefficients.Estimate;
yfit = linspace(0,750,20).*slope;
Rsqrd = mdl.Rsquared.Ordinary;
plot(linspace(0,750,20),yfit,'k-')
ylim([0 750])
xlim([0 750])
xlabel('max fluorescence (AU)')
ylabel('accumulated fluorescence (AU)')
title(['R^2 = ' num2str(Rsqrd)])
%set(gca,'YScale','log'); set(gca,'XScale','log')





%% Max fluo (on X) vs Duration (on Y)
figure
hold on
XforFit = [];
YforFit = [];
for k = 1:length(enhancers)

    % fraction on the X axis
    X = dorsalResultsDatabase.meanAllMaxFluoEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    XError = dorsalResultsDatabase.seAllMaxFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
       % max fluo Y axis
    Y = dorsalResultsDatabase.meanalldurationsEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    YError = dorsalResultsDatabase.sealldurationsEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
    XforFit = [XforFit;X];
    YforFit = [YforFit;Y];
    
    errorbar(X,Y,YError,'o','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
    errorbar(X,Y,XError,'o','horizontal','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
   
end

% remove nans
X = XforFit(~isnan(XforFit));
Y = YforFit(~isnan(XforFit));

%fit to linear model
[p,S] = polyfit(X,Y,1);
f = polyval(p,linspace(0,1,20)); 
%plot(linspace(0,1,20),f,'k-')
% get R^2
[R,P,RL,RU] = corrcoef(X,Y);
Rsqrd = R(1,2);

%fit to line going through the origin
% slope = X\Y; % linear regression throug origin
% yfit = linspace(0,1,20)*slope; % calculate fitted line
mdl = fitlm(X,Y,'Intercept',false);
slope = mdl.Coefficients.Estimate;
yfit = linspace(0,750,20).*slope;
Rsqrd = mdl.Rsquared.Ordinary;
plot(linspace(0,750,20),yfit,'k-')
ylim([0 3.7])
xlim([0 600])
xlabel('maximum spot fluorescence (AU)')
ylabel('spot duration (min)')
title(['R^2 = ' num2str(Rsqrd)])
%set(gca,'YScale','log'); set(gca,'XScale','log')




%% Max fluo (on X) vs Time On (on Y)
figure
hold on
XforFit = [];
YforFit = [];
for k = 1:length(enhancers)

    % fraction on the X axis
    X = dorsalResultsDatabase.meanAllMaxFluoEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    XError = dorsalResultsDatabase.seAllMaxFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
       % max fluo Y axis
    Y = dorsalResultsDatabase.meanTurnOnsEmbryo(...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    YError = dorsalResultsDatabase.seTurnOnsEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}));
    
    XforFit = [XforFit;X];
    YforFit = [YforFit;Y];
    
    errorbar(X,Y,YError,'o','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
    errorbar(X,Y,XError,'o','horizontal','Color',colorList{k},'MarkerFaceColor','w','MarkerEdgeColor',colorList{k},'CapSize',0)
   
end

% remove nans
X = XforFit(~isnan(XforFit));
Y = YforFit(~isnan(XforFit));

%fit to linear model
[p,S] = polyfit(X,Y,1);
f = polyval(p,linspace(0,1,20)); 
%plot(linspace(0,1,20),f,'k-')
% get R^2
[R,P,RL,RU] = corrcoef(X,Y);
Rsqrd = R(1,2);

%fit to line going through the origin
% slope = X\Y; % linear regression throug origin
% yfit = linspace(0,1,20)*slope; % calculate fitted line
mdl = fitlm(X,Y,'Intercept',false);
slope = mdl.Coefficients.Estimate;
yfit = linspace(0,750,20).*slope;
Rsqrd = mdl.Rsquared.Ordinary;
plot(linspace(0,750,20),yfit,'k-')
ylim([0 8])
xlim([0 600])
xlabel('maximum spot fluorescence (AU)')
ylabel('turn on time (min since anaphase)')
title(['R^2 = ' num2str(Rsqrd)])
%set(gca,'YScale','log'); set(gca,'XScale','log')

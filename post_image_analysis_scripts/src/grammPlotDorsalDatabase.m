function grammPlotDorsalDatabase()

close all;

% 
% dataTypes = {'1Dg-12_6_2xDl', 'TwiPEv5(7)_2xDl', '1Dg-20(11)_2xDl','1Dg-8D_FFF', '1Dg11_2xDl', '1DgW_2x_Leica',...
%     '1DgW_FFF', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF', '1Dg-5_2xDl',...
%     '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl', '1DgSVW2_2xDl','1DgVW_2xDl'};


dataTypes = {'1Dg-12_6_2xDl', '1Dg-20(11)_2xDl','1Dg-8D_2xDl' '1Dg11_2xDl', '1Dg-5_2xDl',...
};

affinitiesDataTypes = {'1Dg11_FFF', '1Dg11_2xDl',...'
    '1DgW_FFF', '1DgW_2x_Leica',...
    '1DgVW_FFF'};

phaseDataTypes = {'1Dg11_FFF', '1Dg11_2xDl',...
    '1Dg-5_FFF', '1Dg-8D_FFF', '1Dg-12_6_2xDl', '1Dg-20(11)_2xDl',...
    };

%
% dataTypes = {'1Dg-8D_FFF', '1DgW_2x_Leica',...
%     '1DgW_FFF'};


% try
%     [~, resultsFolder] = getDorsalFolders;
%     load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
% catch
%     dorsalResultsDatabase = createDorsalResultsDatabase(dataTypes);
% end
 dorsalResultsDatabase = createDorsalResultsDatabase(dataTypes);

%% fraction for all data sets, both 1x and 2x (different colors, not concatenated).

figure(2);

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.fracFluoEmbryo;
yLowerError =  dorsalResultsDatabase.meanFracFluoEmbryo -  dorsalResultsDatabase.seFracFluoEmbryo;
yUpperError =  dorsalResultsDatabase.meanFracFluoEmbryo +  dorsalResultsDatabase.seFracFluoEmbryo;


g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12,...
    'color', dorsalResultsDatabase.mother);

% Subdivide the data in subplots horizontally by region of origin
g.facet_grid(dorsalResultsDatabase.enhancer, []);
g.facet_wrap(dorsalResultsDatabase.enhancer, 'ncols', 2);
% Plot raw data as points, then add errorbars
g.geom_point();
% g.geom_interval();
% Plot nonlinear fits of the data with associated confidence intervals

%rate, kd, hill, y offset
y_max = max(dorsalResultsDatabase.fracFluoEmbryo(:));
x_max = max(dorsalResultsDatabase.dorsalFluoBins);
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 2; -y_max];
ub = [y_max*2; Inf; 6; y_max*10];

model = dorsalFitFunction('hill');
modelStr = func2str(model);
modelStr = strrep(modelStr, 'd', 'x');
modelStr = strrep(modelStr, 'p(1)', 'a');
modelStr = strrep(modelStr, 'p(2)', 'b');
modelStr = strrep(modelStr, 'p(3)', 'c');
modelStr = strrep(modelStr, 'p(4)', 'd');
modelStr = strrep(modelStr, '@(p,', '@(a,b,c,d,');
model = str2func(modelStr);

stat_fit(g,'fun', model,'disp_fit', false, 'fullrange', true, ...
    'StartPoint', p0, 'Lower', lb, 'Upper', ub, 'intopt', 'functional', 'geom', 'line');

%set axis labels and legend labels
g.set_names('x','Dorsal Concentration (au)','y','Fraction Active', 'row', '', 'column', '');
grammFigurePBoC(g);
g.draw()

%% max fluo, just phase shifted

figure(5)

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanAllMaxFluoEmbryo;
yLowerError =  dorsalResultsDatabase.meanAllMaxFluoEmbryo -  dorsalResultsDatabase.seAllMaxFluoEmbryo;
yUpperError =  dorsalResultsDatabase.meanAllMaxFluoEmbryo +  dorsalResultsDatabase.seAllMaxFluoEmbryo;

% g=gramm('x',x,...
%     'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
%     'subset', dorsalResultsDatabase.nc==12 &...
%     ...
%    ( strcmpi(dorsalResultsDatabase.enhancer, '1Dg11') |...
%    strcmpi(dorsalResultsDatabase.enhancer, '1Dg-5') |...
%    strcmpi(dorsalResultsDatabase.enhancer, '1Dg-8D') ),...
%     ...
%     'color', dorsalResultsDatabase.mother);
g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ...
   ~strcmpi(dorsalResultsDatabase.enhancer, '1DgW') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVW')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgAW3') &...
      ~strcmpi(dorsalResultsDatabase.enhancer, '1DgS2') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVVW3')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgSVW2'),...
    ...
    'color', dorsalResultsDatabase.mother);


% Subdivide the data in subplots horizontally by region of origin
g.facet_grid(dorsalResultsDatabase.enhancer, [])
% g.facet_wrap(dorsalResultsDatabase.enhancer, 'ncols', 2);
% Plot raw data as points, then plot errorbars
g.geom_point()
g.geom_interval();
% Plot nonlinear fits of the data with associated confidence intervals

%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 2; -y_max];
ub = [y_max*2; Inf; 6; y_max*10];

model = standardizeModelForGramm(dorsalFitFunction('hill'));

stat_fit(g,'fun', model,'disp_fit', false, 'fullrange', true, ...
    'StartPoint', p0, 'Lower', lb, 'Upper', ub, 'geom', 'line');


%set axis labels and legend labels
g.set_names('x','Dorsal Concentration (au)','y','Max fluorescence per trace (au)', 'row', '', 'column', '');
grammFigurePBoC(g);
g.draw()

%% fraction active, phase shifted

figure(6)

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanFracFluoEmbryo;
yLowerError =  dorsalResultsDatabase.meanFracFluoEmbryo -  dorsalResultsDatabase.seFracFluoEmbryo;
yUpperError =  dorsalResultsDatabase.meanFracFluoEmbryo +  dorsalResultsDatabase.seFracFluoEmbryo;

g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ...
   ~strcmpi(dorsalResultsDatabase.enhancer, '1DgW') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVW')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgAW3') &...
      ~strcmpi(dorsalResultsDatabase.enhancer, '1DgS2') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVVW3'),...
    ...
    'color', dorsalResultsDatabase.mother);

% Subdivide the data in subplots horizontally by region of origin
g.facet_grid(dorsalResultsDatabase.enhancer, [])
% g.facet_wrap(dorsalResultsDatabase.enhancer, 'ncols', 2);
% Plot raw data as points, then plot errorbars
g.geom_point()
g.geom_interval();
% Plot nonlinear fits of the data with associated confidence intervals

%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
model = standardizeModelForGramm(dorsalFitFunction('hill'));
stat_fit(g,'fun', model,'disp_fit', false, 'fullrange', true, ...
    'StartPoint', [y_max; x_max/2 ; 1 ; 0],...
    'Lower', [0; 1000; 2; -y_max], 'Upper', [y_max*2; Inf; 6; y_max*10], 'geom', 'line');


%set axis labels and legend labels
g.set_names('x','Dorsal Concentration (au)','y','Fraction active', 'row', '', 'column', '');
grammFigurePBoC(g);
g.draw()

%% fraction active, affinities

figure(7)

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanFracFluoEmbryo;
yLowerError =  dorsalResultsDatabase.meanFracFluoEmbryo -  dorsalResultsDatabase.seFracFluoEmbryo;
yUpperError =  dorsalResultsDatabase.meanFracFluoEmbryo +  dorsalResultsDatabase.seFracFluoEmbryo;

g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ~strcmpi(dorsalResultsDatabase.enhancer, '1Dg-5') &...
    ~strcmpi(dorsalResultsDatabase.enhancer, '1Dg-8D'),...
    'color', dorsalResultsDatabase.mother);

% Subdivide the data in subplots horizontally by region of origin
% g.facet_grid(dorsalResultsDatabase.enhancer, [])
% g.facet_wrap(dorsalResultsDatabase.enhancer, 'ncols', 2);
% Plot raw data as points, then plot errorbars
g.geom_point()
g.geom_interval();
% Plot nonlinear fits of the data with associated confidence intervals

%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 2; -y_max];
ub = [y_max*2; Inf; 6; y_max*10];

model = standardizeModelForGramm(dorsalFitFunction('hill'));

stat_fit(g,'fun', model,'disp_fit', false, 'fullrange', true, ...
    'StartPoint', p0, 'Lower', lb, 'Upper', ub, 'geom', 'line');


%set axis labels and legend labels
g.set_names('x','Dorsal Concentration (au)','y','Fraction active', 'row', '', 'column', '');
grammFigurePBoC(g);
g.draw()

%% time on, affinities

figure;

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanTurnOnsEmbryo;
yLowerError =  dorsalResultsDatabase.meanTurnOnsEmbryo -  dorsalResultsDatabase.seTurnOnsEmbryo;
yUpperError =  dorsalResultsDatabase.meanTurnOnsEmbryo +  dorsalResultsDatabase.seTurnOnsEmbryo;

g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ~strcmpi(dorsalResultsDatabase.enhancer, '1Dg-5') &...
    ~strcmpi(dorsalResultsDatabase.enhancer, '1Dg-8D'),...
    'color', dorsalResultsDatabase.mother);

% Subdivide the data in subplots horizontally by region of origin
g.facet_grid(dorsalResultsDatabase.enhancer, [])
% g.facet_wrap(dorsalResultsDatabase.enhancer, 'ncols', 2);
% Plot raw data as points, then plot errorbars
g.geom_point()
g.geom_interval();
% Plot nonlinear fits of the data with associated confidence intervals

%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 2; -y_max];
ub = [y_max*2; Inf; 6; y_max*10];

model = standardizeModelForGramm(dorsalFitFunction('hill'));

stat_fit(g,'fun', model,'disp_fit', false, 'fullrange', true, ...
    'StartPoint', p0, 'Lower', lb, 'Upper', ub, 'geom', 'line');


%set axis labels and legend labels
g.set_names('x','Dorsal Concentration (au)','y','Time on (min)', 'row', '', 'column', '');
g.axe_property('ylim', [2, 6]);
grammFigurePBoC(g);
g.draw()

%% max fluo, affinities

figure;

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanAllMaxFluoEmbryo;
yLowerError =  dorsalResultsDatabase.meanAllMaxFluoEmbryo -  dorsalResultsDatabase.seAllMaxFluoEmbryo;
yUpperError =  dorsalResultsDatabase.meanAllMaxFluoEmbryo +  dorsalResultsDatabase.seAllMaxFluoEmbryo;

g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ~strcmpi(dorsalResultsDatabase.enhancer, '1Dg-5') &...
    ~strcmpi(dorsalResultsDatabase.enhancer, '1Dg-8D'),...
    'color', dorsalResultsDatabase.mother);

% Subdivide the data in subplots horizontally by region of origin
g.facet_grid(dorsalResultsDatabase.enhancer, [])
% g.facet_wrap(dorsalResultsDatabase.enhancer, 'ncols', 2);
% Plot raw data as points, then plot errorbars
g.geom_point()
g.geom_interval();
% Plot nonlinear fits of the data with associated confidence intervals

%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 2; -y_max];
ub = [y_max*2; Inf; 6; y_max*10];

model = standardizeModelForGramm(dorsalFitFunction('hill'));

stat_fit(g,'fun', model,'disp_fit', false, 'fullrange', true, ...
    'StartPoint', p0, 'Lower', lb, 'Upper', ub, 'geom', 'line');


%set axis labels and legend labels
g.set_names('x','Dorsal Concentration (au)','y','Max fluorescence per trace (au)', 'row', '', 'column', '');
grammFigurePBoC(g);
g.draw()

end
dataTypes = {'1Dg-8D_FFF', '1DgSVW2_2xDl', '1Dg11_2xDl', '1DgW_2x_Leica',...
    '1DgW_FFF', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF', '1Dg-5_2xDl',...
    '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl'};

dataTypes_2x_affinity = {'1Dg11_2xDl', '1DgW_2x_Leica',...
    '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1DgSVW_2xDl', '1DgSVW2_2xDl'};

dataTypes_2x_affinity = {'1Dg11_2xDl', '1DgW_2x_Leica',...
    '1DgVW_FFF', ...
    '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1DgSVW_2xDl', '1DgSVW2_2xDl'};

phaseCell = {'1Dg-20(11)','1Dg-8D', '1Dg-5','1Dg11'};

 [~, resultsFolder] = getDorsalFolders;
    load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
    
close all;    
    %% 
figure(1);

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
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVVW3')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgSVW2')...
    & ~strcmpi(dorsalResultsDatabase.mother, 'FFF'));


% Subdivide the data in subplots horizontally by region of origin
g.facet_grid([], dorsalResultsDatabase.enhancer)

g.set_order_options('column', phaseCell)

g.geom_line()
g.set_names('x','Dorsal Concentration (au)', 'y', 'fraction active',  'row', '', 'column', '');

%default error bars are shaded. this makes them lines
g.geom_interval('geom','errorbar','dodge',0.2,'width',0.8);


%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
model = standardizeModelForGramm(dorsalFitFunction('hill'));
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 1; -y_max];
ub = [y_max*2; Inf; 6; y_max*10];
% stat_fit(g,'fun', model,'disp_fit', true, 'fullrange', true, ...
%     'StartPoint', [y_max; x_max/2 ; 1 ; 0],...
%     'Lower', [0; 1000; 1; -y_max], 'Upper', [y_max*2; Inf; 6; y_max*10], 'geom', 'line');
% stat_fit(g,'fun', model,'disp_fit', true, 'fullrange', true, ...
%     'StartPoint', [y_max; x_max/2 ; 1 ; 0],...
%     'Lower', [0; 1000; 1; -y_max], 'Upper', [y_max*2; Inf; 6; y_max*10]);


grammFigurePBoC(g );
g.draw()

fig = gcf;
ax = fig.Children(2);
ylim(ax, [0, 1]);
xlim(ax, [0, 3250]);


box(ax, 'on')


% standardizeFigure(ax, []);
% ax2 = fig.Children(1);
% standardizeFigure(ax2, []);

%%
figure(2);

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanAllMaxFluoEmbryo;
yLowerError =  dorsalResultsDatabase.meanAllMaxFluoEmbryo -  dorsalResultsDatabase.seAllMaxFluoEmbryo;
yUpperError =  dorsalResultsDatabase.meanAllMaxFluoEmbryo +  dorsalResultsDatabase.seAllMaxFluoEmbryo;

g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ...
   ~strcmpi(dorsalResultsDatabase.enhancer, '1DgW') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVW')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgAW3') &...
      ~strcmpi(dorsalResultsDatabase.enhancer, '1DgS2') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVVW3')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgSVW2')...
    & ~strcmpi(dorsalResultsDatabase.mother, 'FFF'));


% Subdivide the data in subplots horizontally by region of origin
g.facet_grid([], dorsalResultsDatabase.enhancer)

g.set_order_options('column', phaseCell)


g.geom_line()
g.set_names('x','Dorsal Concentration (au)', 'y', 'max fluorescence (au)',  'row', '', 'column', '');

%default error bars are shaded. this makes them lines
g.geom_interval('geom','errorbar','dodge',0.2,'width',0.8);


%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 2; -y_max];
ub = [y_max*2; Inf; 6; y_max*10];

model = standardizeModelForGramm(dorsalFitFunction('hill'));

% stat_fit(g,'fun', model,'disp_fit', true, 'fullrange', true, ...
%     'StartPoint', p0, 'Lower', lb, 'Upper', ub, 'geom', 'line');

stat_fit(g,'fun', model,'disp_fit', true, 'fullrange', true, ...
    'StartPoint', p0, 'Lower', lb, 'Upper', ub);

grammFigurePBoC(g);
g.draw()

fig = gcf;
ax = fig.Children(2);
xlim(ax, [0, 3250]);
ylim(ax, [0, 800]);
box(ax, 'on')

%%


%% Accumulated fluorescence

figure(9);

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanallmrnasEmbryo;
yLowerError =  dorsalResultsDatabase.meanallmrnasEmbryo -  dorsalResultsDatabase.seallmrnasEmbryo;
yUpperError =  dorsalResultsDatabase.meanallmrnasEmbryo +  dorsalResultsDatabase.seallmrnasEmbryo;

g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ...
   ~strcmpi(dorsalResultsDatabase.enhancer, '1DgW') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVW')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgAW3') &...
      ~strcmpi(dorsalResultsDatabase.enhancer, '1DgS2') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVVW3')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgSVW2')...
    & ~strcmpi(dorsalResultsDatabase.mother, 'FFF'));


% Subdivide the data in subplots horizontally by region of origin
g.facet_grid([], dorsalResultsDatabase.enhancer)

g.set_order_options('column', phaseCell)


g.geom_line()
g.set_names('x','Dorsal Concentration (au)', 'y', 'accumulated fluorescence (au)',  'row', '', 'column', '');

%default error bars are shaded. this makes them lines
g.geom_interval('geom','errorbar','dodge',0.2,'width',0.8);


%rate, kd, hill, y offset
y_max = nanmax(y(:));
x_max = max(x);
p0 = [y_max; x_max/2 ; 1 ; 0];
lb = [0; 1000; 2; -y_max];
ub = [Inf; Inf; 6; y_max*10];


model = standardizeModelForGramm(dorsalFitFunction('hill'));


% stat_fit(g,'fun', model,'disp_fit', true, 'fullrange', true, ...
%     'StartPoint', p0, 'Lower', lb, 'Upper', ub, 'geom', 'line');

stat_fit(g,'fun', model,'disp_fit', true, 'fullrange', true, ...
    'StartPoint', p0, 'Lower', lb, 'Upper', ub);

grammFigurePBoC(g);
g.draw()

fig = gcf;
ax = fig.Children(2);
xlim(ax, [0, 3250]);
ylim(ax, [0, 1250]);
box(ax, 'on')
%%
%%

figure(10);

x = dorsalResultsDatabase.dorsalFluoBins;
y = dorsalResultsDatabase.meanTurnOnsEmbryo;
yLowerError =  dorsalResultsDatabase.meanTurnOnsEmbryo -  dorsalResultsDatabase.seTurnOnsEmbryo;
yUpperError =  dorsalResultsDatabase.meanTurnOnsEmbryo +  dorsalResultsDatabase.seTurnOnsEmbryo;

g=gramm('x',x,...
    'y', y, 'ymin', yLowerError , 'ymax', yUpperError,...
    'subset', dorsalResultsDatabase.nc==12 &...
    ...
   ~strcmpi(dorsalResultsDatabase.enhancer, '1DgW') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVW')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgAW3') &...
      ~strcmpi(dorsalResultsDatabase.enhancer, '1DgS2') &...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgVVW3')&...
     ~strcmpi(dorsalResultsDatabase.enhancer, '1DgSVW2')...
    & ~strcmpi(dorsalResultsDatabase.mother, 'FFF'));


% Subdivide the data in subplots horizontally by region of origin
g.facet_grid([], dorsalResultsDatabase.enhancer)

g.set_order_options('column',phaseCell)


g.geom_line()
g.set_names('x','Dorsal Concentration (au)', 'y', 'activation time (min)',  'row', '', 'column', '');

%default error bars are shaded. this makes them lines
g.geom_interval('geom','errorbar','dodge',0.2,'width',0.8);

% grammFigurePBoC(g);
g.draw()

fig = gcf;
ax = fig.Children(2);
xlim(ax, [0, 3250]);
% ylim(ax, [0, 1250]);
box(ax, 'on')
%%
% 

scores = [-20,-8, -5, 0]';

dorsalResultsDatabase.patserScore = nan(length(dorsalResultsDatabase.enhancer), 1);
for k = 1:length(phaseCell)
    for j = 1:length(dorsalResultsDatabase.enhancer)
        if strcmpi(dorsalResultsDatabase.enhancer{j}, phaseCell{k})
            dorsalResultsDatabase.patserScore(j) = scores(k);
        end
    end
end


maxFraction = nan(length(scores), 1);
meanFraction = nan(length(scores), 1);
maxFluo = nan(length(scores), 1);
meanFluo = nan(length(scores), 1);
sumFluo = nan(length(scores), 1);
sumFluoSE = nan(length(scores), 1);
sumFluoFracProduct = nan(length(scores), 1);
sumFluoFracProductSE = nan(length(scores), 1);

for k = 1:length(scores)
    
    relation_frac = dorsalResultsDatabase.meanFracFluoEmbryo( dorsalResultsDatabase.patserScore == scores(k)  &...
       strcmpi(dorsalResultsDatabase.mother,'2x') );
   
   relation_fracSE = dorsalResultsDatabase.seFracFluoEmbryo( dorsalResultsDatabase.patserScore == scores(k)  &...
       strcmpi(dorsalResultsDatabase.mother,'2x') );
   
    relation_fluo = dorsalResultsDatabase.meanallmrnasEmbryo( dorsalResultsDatabase.patserScore == scores(k)  &...
       strcmpi(dorsalResultsDatabase.mother,'2x') );
   
   relation_fluoSE = dorsalResultsDatabase.seallmrnasEmbryo( dorsalResultsDatabase.patserScore == scores(k)  &...
       strcmpi(dorsalResultsDatabase.mother,'2x') );
    
   try
        maxFraction(k) = nanmax(relation_frac);
        meanFraction(k) = nanmean(relation_frac); 
        maxFluo(k) = nanmax(relation_fluo);
        meanFluo(k) = nanmean(relation_fluo);
        
        sumFluo(k) = nansum(relation_fluo);
        sumFluoSE(k) = sqrt( nansum(relation_fluoSE.^2 ) );
        
        sumFluoFracProduct(k) = nansum(relation_fluo.*relation_frac);
        %crossing my fingers this next line doesn't have mistakes
        sumFluoFracProductSE(k) = sqrt( nansum(...
            (relation_fluo.*relation_frac) .* sqrt( (relation_fluoSE./relation_fluo).^2 + (relation_fracSE./relation_frac).^2 )...
            ).^2 );


   catch
        %fails when results are empty for that dataset. ie it hasn't been
        %pushed through the pipeline.
        maxFraction(k) = nan;
        meanFraction(k) = nan;
        maxFluo(k) = nan;
        meanFluo(k) = nan;
        sumFluo(k) = nan;
        sumFluoSE(k) = nan;
        sumFluoFracProduct(k) = nan;
   end
     
end

figure(3)
tiledlayout('flow');

nexttile
plot(scores, maxFraction);
xlim([4.5, 6.5])
ylim([0, 1])
xlabel('binding site affinity (patser score)')
ylabel('max fraction active')
standardizeFigure(gca, []);

nexttile
plot(scores, meanFraction);
xlim([4.5, 6.5])
ylim([0, 1])
xlabel('binding site affinity (patser score)')
ylabel('mean fraction active')
standardizeFigure(gca, []);

nexttile
plot(scores, maxFluo);
xlim([4.5, 6.5])
% ylim([0, 1])
xlabel('binding site affinity (patser score)')
ylabel('max fluorescence (au)')
standardizeFigure(gca, []);

nexttile
plot(scores, meanFluo);
xlim([4.5, 6.5])
% ylim([0, 1])
xlabel('binding site affinity (patser score)')
ylabel('mean fluorescence (au)')
standardizeFigure(gca, []);

nexttile
errorbar(scores, sumFluo, sumFluoSE);
xlim([4.5, 6.5])
% ylim([0, 1])
xlabel('binding site affinity (patser score)')
ylabel('accumulated mean fluorescence (au)')
standardizeFigure(gca, []);
title('Integrated across DV. Mean across time. Taken only across active nuclei')


nexttile
errorbar(scores, sumFluoFracProduct, sumFluoFracProductSE);
xlim([4.5, 6.5])
% ylim([0, 1])
xlabel('binding site affinity (patser score)')
ylabel('accumulated mean fluorescence (au)')
standardizeFigure(gca, []);
title(gca, 'Integrated across DV. Mean across time. Taken across active AND inactive nuclei.')


%%
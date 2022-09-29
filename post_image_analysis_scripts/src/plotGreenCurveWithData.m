function [dat_onset, dat_fraction] = plotGreenCurveWithData(varargin)

if nargout > 0
    displayFigures = false;
else
    displayFigures = true;
end
    
%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

numBins = 20;

fiducialTime = 6; %mins?

for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

%just gonna grab 1dg11_2xDl data
enhancerStruct = combinedCompiledProjects_allEnhancers(...
    [combinedCompiledProjects_allEnhancers.cycle]==12 &...
    ({combinedCompiledProjects_allEnhancers.dataSet} == "1Dg11_2xDl" |...
      {combinedCompiledProjects_allEnhancers.dataSet} == "1Dg11_FFF"  )&...
    ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));


% bin nuclei
%this function calculates the Dorsal fluorescence at some arbitrary time
%in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);
nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];

binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end

coveredBins = unique([enhancerStruct.dorsalFluoBin2]);
binValues = binValues(coveredBins);
%store everything in these arrays
mean_fraction_acrossNuclei_perBin = [];
mean_fraction_acrossEmbryos_perBin = [];
se_fraction_acrossEmbryos_perBin = [];

OnNucleiPerBin = [];
OffNucleiPerBin = [];

prefixes = unique({enhancerStruct.prefix});
dat_fraction = [];
dat_onset = [];
for p = 1:length(prefixes)
    
    mean_fraction_acrossNuclei_perBin = [];
    mean_timeOn_acrossNuclei_perBin = [];
    se_timeOn_acrossNuclei_perBin = [];
    
    for b = 1:length(coveredBins)
        
        binID = coveredBins(b);
        binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin2]== binID & strcmpi({enhancerStruct.prefix}, prefixes{p}) );
        
        activeNuc_Bin = length([binStruct.particleTimeOn]);
        inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;
        mean_fraction_acrossNuclei_perBin(b) = activeNuc_Bin/length(binStruct);
        mean_timeOn_acrossNuclei_perBin(b) = nanmean([binStruct.particleTimeOn]);
        se_timeOn_acrossNuclei_perBin(b) = nanstd([binStruct.particleTimeOn])./sqrt(activeNuc_Bin);
    end
    dat_fraction = [dat_fraction, mean_fraction_acrossNuclei_perBin];
    dat_onset = [dat_onset, mean_timeOn_acrossNuclei_perBin];
    
    dat_fraction_dl(p, :) = mean_fraction_acrossNuclei_perBin;
    dat_onset_dl(p, :) = mean_timeOn_acrossNuclei_perBin;
    dat_onsetSE_dl(p, :) = mean(se_timeOn_acrossNuclei_perBin);
    
end

% save([resultsFolder, filesep, '2Dhist.mat'], 'dat_onset', 'dat_fraction')

 dat_fraction_dl0 =  dat_fraction_dl;
 dat_onset_dl0 = dat_onset_dl;

%clean data
dat_fraction_dl(dat_onset_dl<2) = nan;
dat_onset_dl(dat_onset_dl < 2) = nan;
dat_fraction_dl(dat_onset_dl>8) = nan;
dat_onset_dl(dat_onset_dl>8) = nan;

dat_fraction_dl_mean = nanmean(dat_fraction_dl,1);
dat_fraction_dl_ste = nanstd(dat_fraction_dl,1)./sqrt(length(prefixes));

dat_onset_dl_mean = nanmean(dat_onset_dl, 1);
dat_onset_dl_ste = nanstd(dat_onset_dl, 1)./sqrt(length(prefixes));

c = binValues; 

dat_onset_dl_mean(isnan(dat_onset_dl_ste)) = []; 
c(isnan(dat_onset_dl_ste)) = [];
dat_fraction_dl_mean(isnan(dat_onset_dl_ste)) = []; 
dat_fraction_dl_ste(isnan(dat_onset_dl_ste)) = []; 
dat_onset_dl_ste(isnan(dat_onset_dl_ste)) = []; 



%this is for scatter_ellipse. not a huge fan, but leaving this here for
%future reference. 
for k = 1:length(c)
    cov(:, :, k) = diag([dat_fraction_dl_ste(k), dat_onset_dl_ste(k)]);
end
 
P = [dat_fraction_dl_mean;dat_onset_dl_mean]';

Ps = sortrows(P, 1);

x = Ps(:, 1);
y = Ps(:, 2);
figure; 
colormap(brewermap(length(c),'Greens'))
scatter(y, x,[],c, 'o', 'filled', 'MarkerEdgeColor', 'k')
hold on
errorbar(y, x, dat_fraction_dl_ste, dat_fraction_dl_ste,...
    dat_onset_dl_ste, dat_onset_dl_ste,...
    'k', 'LineStyle', 'none', 'CapSize', 0)
% h = scatter_ellipse(x, y, c, cov)
ylim([0, 1])
xlim([0, 8.5])
xlabel('mean transcriptional onset time (min)')
ylabel('fraction of active nuclei')
title('1Dg11_2xDl')
colorbar;

%clean data
dat_fraction(dat_onset<2) = [];
dat_onset(dat_onset < 2) = [];
dat_fraction(dat_onset>8) = [];
dat_onset(dat_onset>8) = [];
dat_fraction(isnan(dat_onset)) = [];
dat_onset(isnan(dat_onset)) = [];
dat_onset(isnan(dat_fraction)) = [];
dat_fraction(isnan(dat_fraction)) = [];

P = [dat_onset;dat_fraction]';

Ps = sortrows(P, 2);

x = Ps(:, 2);
y = Ps(:, 1);

figure; plot(x, y,'ok')




nBins = [15, 5];

if displayFigures
    
    fighis = figure;
    ax = axes(fighis);
    
    
    figure;
    colormap(brewermap(20,'Purples'))
%     yyaxis left
%     h = binscatter(dat_onset, dat_fraction, nBins);
%     h.ShowEmptyBins = 'on';
    g = histogram2(dat_onset, dat_fraction,nBins,...
        'DisplayStyle','tile','ShowEmptyBins','on');
    g.EdgeColor = 'none';
    xlim([0, 10])
    ylim([0, 1])
    grid off
    colorbar
    xlabel('mean transcription onset time (min)')
    ylabel('fraction active')
    hold on
%     [k,	~] = convhull([x;y]');
    % plot(x(k),y(k))
%     hold on
    % fill(x(k),y(k), 'g', 'FaceAlpha', .3, 'LineWidth', 3, 'EdgeColor', 'g')
    % hull = alphaShape(x_hull, y_hull,Inf,'HoleThreshold',1E30 );
    % colormap brewermap(20,'Blues')
    g = histogram2(ax, dat_onset, dat_fraction,nBins, 'DisplayStyle','tile','ShowEmptyBins','on',...
        'Normalization', 'probability');
    

    %Now, we'll add the marginalized turn on density
    bincenters = g.XBinEdges + g.BinWidth(1)/2;
    bincenters = bincenters(1:end-1);
    
    %let's integrate out the fraction axis 
    f = sum(g.Values, 2);
    
    hold on
    yyaxis right
    plot(bincenters, f, '-ok', 'LineWidth', 2);
    ylabel('density')
    
    
    %we're going to plot the narrowest region containing 95% density. 
    width = nan(length(f), length(f));
    target = .95;
    for k = 1:length(f)-1
       for j = k+1:length(f)
           ss = sum(f(k:j));
           if ss >= target && j > k
            width(k, j) = abs(j-k);
           end
       end
    end
    [~,I] = min(width, [], 'all', 'linear');
    [row,col] = ind2sub([length(f), length(f)], I);

    xline(bincenters(row), 'g', 'LineWidth', 2)
    hold on
    xline(bincenters(col), 'g', 'LineWidth', 2)
    
    
   
    
%     scatter(x, y)
%     
%%    
    figure;
    
    


end

%%



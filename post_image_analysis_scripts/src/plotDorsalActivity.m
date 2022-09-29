function [fit, model] = plotDorsalActivity(x, y,activity, nc,...
    DataType, ymean, se, plotScatter, varargin)

    arguments
       x double
       y double
       activity string
       nc double
       DataType string
       ymean double
       se double
       plotScatter logical      
    end
    
    arguments(Repeating)
        varargin
    end
    fitOpts = varargin;

    xMesh = repmat(x, size(y,2), 1)';
    [fit, model] = fitDorsalActivity(xMesh, y, fitOpts{:});

    idx = ~any(isnan(ymean),2); %remove nans from x mesh
    xFiltered = xMesh(idx);
    plotRange = min(xFiltered(:)) : 1: xFiltered(end)*1.3;
    
    figure('Units', 'points', 'Position', [0, 0, 200, 200]);
    clr = 'r';
    plot(plotRange, model(fit,plotRange),...
        '-', 'DisplayName',['fit: ',num2str(round(fit))],...
        'LineWidth', 2, 'Color', clr);
    set(gca,'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
    axis square
    hold on
    
    if plotScatter
     errorbar(xFiltered, ymean(idx), se(idx), 'o', 'DisplayName',...
         DataType, 'MarkerSize', 4, 'MarkerFaceColor', clr, 'CapSize', 0);
    end
    
    xlabel('dorsal concentration (au)');
    ylabel(activity);
    
    
    if nc < 6
        nc = nc + 11;
    end
    title({DataType, [' nc: ', num2str(nc)]});

    legend(DataType, { num2str(round(fit)), 'amplitude, KD(au), n, y offset'} )

    leg = get(gca, 'Legend'); w=.02;h=.01;
    set(leg, 'Units', 'normalized', 'Position', [0.5,0.8,w, h], 'Box','off');

end
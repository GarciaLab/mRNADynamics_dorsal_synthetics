
lims = [0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995];
figure;
tiledlayout(1, numel(enhancers));

clr = viridis(length(out.predlims));
y_lower = {};
y_upper = {};
y_median = {};
yl = [];
yu = [];
ym = [];
nFVals = {};

for k = 1:length(out.predlims) %%iterate over all enhancers
    nn = (size(out.predlims{k}{1},1) + 1) / 2;
    plimi = out.predlims{k};
    %     x{k} = data{k}.ydata(:, 1);
    
    nexttile;
    for j = 1:length(plimi)
        
        y_lower{k, j} = plimi{j}(1,:);
        %         y_lower{k, j} = plimi{j}(4,:);
        y_upper{k, j}  = plimi{j}(2*nn-1,:);
        %         y_upper{k, j}  = plimi{j}(6,:);
        
        y_median{k, j} = plimi{j}(nn,:);
        
%         yl{k, j} = nan(1, length(xForPlot));
%         yu{k, j} = nan(1, length(xForPlot));
%         ym{k, j} = nan(1, length(xForPlot));
        
        yl{k, j} = y_lower{k, j};
        yu{k, j} = y_upper{k, j};
        ym{k, j} = y_median{k, j};
%         for i = 1:length(xForPlot)
%             yl{k, j}(i) = nanmedian(y_lower{k, j});
%             yu{k, j}(i) = nanmedian(y_upper{k, j});
%             ym{k, j}(i) = nanmedian(y_median{k, j});
%         end
        for i = 1:length(binMidValues)
            temp_f = MaxFluosPerEmbryoAll{k}(:, i);
            temp_f(temp_f == 0) = nan;
            nFVals{k}(i) = sqrt(length(temp_f(~isnan(temp_f))));
        end
        
        %Plot the prediction intervals
        yl2 =  yl{k, j};
        yu2 = yu{k, j};
        yl2(isnan(yl2)) = [];
        yu2(isnan(yu2)) = [];
        fillyy(xForPlot, yl2, yu2,[0.9 0.9 0.9]);
        hold on
        
        %Plot the data with errorbars
        nFVals{k}(nFVals{k}==0) = nan;
        yyy = MaxFluosPerEmbryoAll{k};
        yyy(yyy==0) = nan;
        y4 = nanmean(yyy, 1);
        ye4 =  nanstd(yyy, 1)./nFVals{k};
        x4 = binMidValues(~isnan(y4));
        y4(isnan(y4)) = [];
        ye4(isnan(ye4)) = [];
        errorbar(x4, y4, ye4,...
            'o', 'CapSize', 0, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', clr(k, :) )
        ylabel('max. fluo.(a.u.)')
        xlabel('Dorsal concentration (a.u.)')
        ylim([50, 450]);
        xlim([80, 1E5]);
        
        %Plot the prediction curve
        
        %Plot only over data
        %         y5 = ym{k, j};
        %         x5 = binMidValues(~isnan(y5));
        %         y5 = y5(~isnan(y5));
        %         plot(x5, y5, '-')
        
        %Plot only over data
        y5 = ym{k, j};
        %             x5 = binMidValues(~isnan(y5));
        %             y5 = y5(~isnan(y5));
        plot(xForPlot, y5, '-')
        
        %Extrapolate past data
        %         plot(xForPlot, yy{k}, '-', 'LineWidth', 2)
        set(gca, 'XScale', 'log')
        set(gca, 'XTick',[0, 100, 1E3, 1E4, 1E5])
        
        
%         axis square;
    end
    
end
function t = plotModelCoeffs(fits, scores, enhancers)

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

xlabel(t, 'affinity (Patser score)')
ylabel(t, 'predicted value')

end
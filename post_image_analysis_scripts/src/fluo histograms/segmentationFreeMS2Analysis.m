function out = segmentationFreeMS2Analysis(Prefix, nBins, varargin)

%color is the bin color accepted as a string within one element cell.
%use an empty cell as default. accepted colors- 'red', 'yellow', 'cyan', 'magenta',
%'lightBlue'

scale = 1;
color = {'red'}; %just using red as the default
optionalResults = '';
ax = [];
shouldMaskNuclei = false;
vals = [];
out = struct;

%area for improvement: segment the nuclei and only include dog values from
%within nuclei 
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'color')
        color = varargin{i+1};
    elseif strcmpi(varargin{i}, 'scale')
        scale = varargin{i+1};
    elseif strcmpi(varargin{i}, 'ax')
        ax = varargin{i+1};
    elseif strcmpi(varargin{i}, 'nuclearMask')
        shouldMaskNuclei = true;
    elseif strcmpi(varargin{i}, 'vals')
        vals = varargin{i+1};
    end
end

if isempty(vals)
    vals = getDoGVals(Prefix, shouldMaskNuclei);
end


if isempty(ax)
    fig = figure(); %#ok<NASGU>
    ax = axes(figure);
end

h = histogram(ax, vals,nBins,'Normalization','pdf', 'facealpha', .6);
hold on;
pd = fitdist(vals','Normal');
x_values = 0:.001:.2;
y = pdf(pd,x_values);
plot(ax, x_values,y,'-','LineWidth',.5)
% h = histfit(ax, vals,nBins);
% edges = 0:.05:.4;
% h = histogram(ax, vals, edges,'Normalization','pdf', 'facealpha', .6);
% set(ax,'YScale','log');
xlabel(ax, 'log(max DoG intensity + 1) (au)');
ylabel(ax, 'frequency');
title(Prefix, 'Interpreter', 'none');
% standardizeFigure(ax, [], color{1}, 'fontSize', 14);

out.xvalues = x_values;
out.y = y;
out.vals = vals;
out.pd = pd;
out.h = h;
out.ax = ax;

end

% load('C:\Users\Armando\Dropbox\DorsalSyntheticsDropbox\tfentryexit_paramsearch.mat')
function goodMatrixIndices = plotTFDrivenParams_3params(factive, mfpts, varargin)

nPoints = []; %if this is an integer, plotting will subsample using this many points. Otherwise, no subsampling
fig = [];
goodMatrixIndices = [];
params = struct;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


mfpts0 = mfpts;
factive0 = factive;

factive(isnan(mfpts)) = [];
mfpts(isnan(mfpts)) = [];

x = factive(:);
y = mfpts(:);

%% make  convex hull out of the allowable region
% fmin = .1;
fmin = 0;
fmax = 1;
tmin = 3.1; %95% highest density interval
tmax = 7.1; %95% highest density interval
deltatmin = 0; %arbitrary
deltatmax = 3; %arbitrary

box = makeBox([fmin, tmin, deltatmin], [fmax, tmax, deltatmax]);
x_hull = box(:, 1);
y_hull = box(:, 2);



hull_2D = alphaShape(x_hull, y_hull, Inf,'HoleThreshold',1E30 );
%%

%% Subsample and/or round
if ~isempty(nPoints)
    if length(x) > nPoints
        subSample = randsample(length(x),nPoints);
        x = x(subSample);
        y = y(subSample);
    end
end


in_2D = inShape(hull_2D,x, y);
in = in_2D;

%%
if isempty(fig)
    fig = figure();
end


if nargout == 0
    
    dim_dl = size(factive0, 1);
    
    factive0(isnan(mfpts)) = nan;
    mfpts0(isnan(mfpts)) = nan;
    
    %these points are unreliable due to small number issues,
    %so let's remove them.
    mfpts0(factive0 < .05) = nan;
    factive0(factive0 < .05) = nan;
    
    x0 = x;
    y0 = y;
    x0(x0 < .05) = nan;
    y0(x0 < .05) = nan;
    
    %%
    
    scatter(y0,x0,'o', 'MarkerFaceColor', [128 128 128]/255,...
        'MarkerEdgeColor', 'none')
    hold on
    
    scatter(y0(in),x0(in),'o', 'MarkerFaceColor', 'r',...
        'MarkerEdgeColor', 'none')
    hold on
    
    colormap(brewermap(dim_dl,'Greens'))
  try  
    i = nearestIndex(params.pi_basics, 2);
    j = nearestIndex(params.pi_exits, .001);
    k = nearestIndex(params.pi_entries, 2);
  end
    
%     f = @(t) reshape( t(i, j, k), [numel( t(i, j, k) ), 1]);
    
    
%            scatter( f(factive0), f(mfpts0), [],params.dls, 'o', 'filled');
    
%     
%     yy = f(factive0)';
%     xx = f(mfpts0)';
%     zz = zeros(size(xx));
%     col = params.dls; % This is the color, vary with x in this case.
%     surface([xx;xx],[yy;yy],[zz;zz],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',3);
    
%     colorbar;
    ylim([0, 1])
    xlim([0, 10])
    %             set(gca,'Color','r')
    
    hold on
    %%
    
    set(gca, 'ColorScale', 'linear');
    ylim([0, 1]);
    xlim([0, 8.5]);
%     ylabel('fraction of active nuclei')
%     xlabel('mean transcription onset time (min)')
    
    
    %         set(topAxes,'xtick',[],'ytick',[]);
    
    pbaspect([3 1 1])
    
else
    %%% Let's return the good parameters
    if   isempty(nPoints)
        
        %let's change nans to some numerical values that will get rejected
        factive0(isnan(mfpts0)) = -1;
        mfpts0(isnan(mfpts0)) = 100;
        
        
        
        in_temp_2D = inShape(hull_2D, factive0(:), mfpts0(:) );
        goodLinearIndices = find(in_temp_2D);
        [in1, in2, in3] = ind2sub(size(mfpts0), goodLinearIndices);
        goodMatrixIndices = [in1 in2 in3];
        
        
        
        %dls, kds, pi1s, cs, pi2s
        fields_params = fieldnames(params);
        figure
        
        t = tiledlayout('flow');
        title(t, 'x axes all log10 scale')
        fields_params = {'pi_basics', 'pi_entries', 'pi_exits'};
        for k = 1:size(goodMatrixIndices, 2)
            nexttile;
%             try
                histogram(log10(params.(fields_params{k})(goodMatrixIndices(:, k))));
%             end
            xlabel(fields_params{k});
            
        end
    end
    
end

end
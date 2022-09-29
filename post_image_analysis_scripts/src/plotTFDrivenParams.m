% load('C:\Users\Armando\Dropbox\DorsalSyntheticsDropbox\tfentryexit_paramsearch.mat')
function goodMatrixIndices = plotTFDrivenParams(factive, dt, mfpts, varargin)

nPoints = []; %if this is an integer, plotting will subsample using this many points. Otherwise, no subsampling
shouldRound = false; %if true, subsample by discretizing data points and rounding nearby values
fig = [];
goodMatrixIndices = [];
dim = 2;
params = struct;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

[~, resultsFolder] = getDorsalFolders;

mfpts0 = mfpts;
dt0 = dt;
factive0 = factive;

factive(isnan(mfpts)) = [];
dt(isnan(mfpts)) = [];
mfpts(isnan(mfpts)) = [];

factive(isnan(dt)) = [];
mfpts(isnan(dt)) = [];
dt(isnan(dt)) = [];

x = factive(:);
y = mfpts(:);
z = dt(:);

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
z_hull = box(:, 3);


hull = alphaShape(x_hull, y_hull, z_hull,Inf,'HoleThreshold',1E30 );

hull_2D = alphaShape(x_hull, y_hull, Inf,'HoleThreshold',1E30 );
%%

%% Subsample and/or round
if ~isempty(nPoints)
    if length(x) > nPoints
        subSample = randsample(length(x),nPoints);
        x = x(subSample);
        y = y(subSample);
        z = z(subSample);
    end
end

if shouldRound
    r_rounded = unique([round(x,2) round(y,1) round(z,-1)], 'rows'); %factive, mfpts, dt
    x = r_rounded(:, 1);
    y = r_rounded(:, 2);
    z = r_rounded(:, 3);
end

in = inShape(hull,x, y, z);

in_2D = inShape(hull_2D,x, y);

%%
if ~isempty(fig)
    gcf;
end


if nargout == 0
    
    if dim == 3
        plot(hull)
        hold on
        
        scatter3(x(in),y(in),z(in),'r.')
        scatter3(x(~in),y(~in), z(~in),'b.')
        %
        % xlabel('factive')
        % ylabel('mean turn on (min)')
        % zlabel('delta t')
        % legend('viable region', 'viable parameters', 'unphysical parameters');
        
        ax = gca;
        % ax.Children(3).EdgeColor = 'none';
        
        % title('Parameter space for TF Driven model')
        
        axis square;
        
        xlim([0, 1]);
        ylim([0, 10]);
        fig = gcf;
        fig.Renderer='Painters';
        
        view(0, -90)
        
    elseif dim == 2
        
        %         bottomAx =axes;
        
        %         colormap(brewermap(20,'Blues'));
        
        %         try
        %             load([resultsFolder, filesep, '2Dhist.mat'], 'dat_onset', 'dat_fraction')
        %         catch
        %             [dat_onset, dat_fraction] = plotGreenBoxWithData;
        %         end
        
        %         nBins = [6, 16];
        %         h = binscatter(dat_fraction,dat_onset, nBins);
        %         h.ShowEmptyBins = 'on';
        %         colormap(brewermap(20,'Blues'));
        xlim([0, 1])
        ylim([0, 10]);
        
        xlabel('fraction of active nuclei')
        ylabel('mean transcription onset time (min)')
        %          legend('viable region', 'viable parameters', 'unphysical parameters');
        
        hold on
        
        %         scatter(x(in),y(in),'r.')
        %         scatter(x(in_2D),y(in_2D),'r.')
        %         scatter(x(in_2D),y(in_2D),'k.')
        
        hold on
        
        %         scatter(x(~in),y(~in),'b.')
        %         scatter(x(~in_2D),y(~in_2D),'b.')
        %         scatter(x(~in_2D),y(~in_2D),'k.')
        
        %         set (gca,'Ydir','reverse')
        
        %          [in1, in2, in3, in4, in5] = ind2sub(size(mfpts0), goodLinearIndices)
        dim_dl = size(factive0, 1);
        %         c = brewermap(dim_dl,'Reds');
        %         c = brewermap(dim_dl,'Reds');
        
        
        factive0(isnan(mfpts)) = [];
        dt0(isnan(mfpts)) = [];
        mfpts0(isnan(mfpts)) = [];
        
        %these points are unreliable due to small number issues,
        %so let's remove them.
        mfpts0(factive0 < .05) = nan;
        dt0(factive0 < .05) = nan;
        factive0(factive0 < .05) = nan;
        
        %dls, kds, pi1s, cs, pi2s
        
        %         figure
        if ~isempty(params)
            j = nearestIndex(params.kds, 1E5); %10. %kd==100,000
            l = 1;
            m = nearestIndex(params.cs, 77);%5; %c == 77
        else
            j=10; %kd==100,000
            l=1;
            m=5;%c == 77
        end
        
        x0 = x;
        y0 = y;
        x0(x0 < .05) = nan;
        y0(x0 < .05) = nan;
        %         scatter(x0(in),y0(in),'o', 'MarkerFaceColor', [128 128 128]/255,...
        %             'MarkerEdgeColor', 'none')
        hold on
        colormap(brewermap(dim_dl,'Greens'))
        
        if params.model == "basic"
            %             scatter(reshape(factive0(:, j, l, m), [numel(factive0(:, j, l, m)), 1]),...
            %                 reshape(mfpts0(:, j, l, m), [numel(mfpts0(:, j, l, m)), 1]),[],...
            %                 params.dls, 'o', 'filled')
            
            
            f = @(t) reshape(t(:, j, l, m), [numel(t(:, j, l, m)), 1]);
            
            xx = f(factive0)';
            yy = f(mfpts0)';
            zz = zeros(size(xx));
            col = params.dls; % This is the color, vary with x in this case.
            surface([xx;xx],[yy;yy],[zz;zz],[col;col],...
                'FaceColor','none',...
                'EdgeColor','interp',...
                'LineWidth',3);
            
            colorbar;
            xlim([0, 1])
            ylim([0, 10])
            
            
            hold on
            scatter(x0(in),y0(in),'o', 'MarkerFaceColor', [128 128 128]/255,...
                'MarkerEdgeColor', 'none')
        else
            %%
            figure
            %dls, kds, pi1s, cs, pi2s
            
            
            scatter(x0(in),y0(in),'o', 'MarkerFaceColor', [128 128 128]/255,...
                'MarkerEdgeColor', 'none')
            hold on
            
            colormap(brewermap(dim_dl,'Greens'))
            if ~params.exitOnlyDuringOffStates
                j = nearestIndex(params.kds, 1E4);
                m = nearestIndex(params.cs, 10);
                l = nearestIndex(params.pi1s, .001);
                n = nearestIndex(params.pi2s, 2);
            else
                j = nearestIndex(params.kds, 1E4);
                m = nearestIndex(params.cs, 10);
                l = nearestIndex(params.pi1s, .001);
                n = nearestIndex(params.pi2s, 2.5);
            end
            
            f = @(t) reshape(t(:, j, l, m, n), [numel(t(:, j, l, m, n)), 1]);
            
            
            %            scatter( f(factive0), f(mfpts0), [],params.dls, 'o', 'filled');
            
            
            xx = f(factive0)';
            yy = f(mfpts0)';
            zz = zeros(size(xx));
            col = params.dls; % This is the color, vary with x in this case.
            surface([xx;xx],[yy;yy],[zz;zz],[col;col],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',3);
            
            colorbar;
            xlim([0, 1])
            ylim([0, 8.5])
            %             set(gca,'Color','r')
            
            hold on
            
            hold on
            %%
        end
        
        set(gca, 'ColorScale', 'linear');
        xlim([0, 1]);
        ylim([0, 10]);
        
        %         set(topAxes,'xtick',[],'ytick',[]);
    end
else
    %%% Let's return the good parameters
    if ~shouldRound && isempty(nPoints)
        
        %let's change nans to some numerical values that will get rejected
        factive0(isnan(mfpts0)) = -1;
        dt0(isnan(mfpts0)) = 100;
        mfpts0(isnan(mfpts0)) = 100;
        
        
        if dim==3
            
            factive0(isnan(dt0)) = -1;
            mfpts0(isnan(dt0)) = 100;
            dt0(isnan(dt0)) = 100;
            
            in_temp = inShape(hull,factive0(:), mfpts0(:), dt0(:));
            goodLinearIndices = find(in_temp);
            [in1, in2, in3, in4, in5] = ind2sub(size(mfpts0), goodLinearIndices);
            goodMatrixIndices = [in1 in2 in3 in4 in5];
        elseif dim==2
            
            in_temp_2D = inShape(hull_2D, factive0(:), mfpts0(:) );
            goodLinearIndices = find(in_temp_2D);
            [in1, in2, in3, in4, in5] = ind2sub(size(mfpts0), goodLinearIndices);
            goodMatrixIndices = [in1 in2 in3 in4 in5];
            
        end
        
        
        %dls, kds, pi1s, cs, pi2s
        fields_params = fieldnames(params);
        figure;
        
        t = tiledlayout('flow');
        title(t, 'x axes all log10 scale')
        
        for k = 1:size(goodMatrixIndices, 2)
            nexttile;
            try
                histogram(log10(params.(fields_params{k})(goodMatrixIndices(:, k))));
            end
            xlabel(fields_params{k});
            
        end
    end
    
end

end
function h = scatter_ellipse(X,Y,C,Cov,varargin)
%SCATTER_ELLIPSE - colored ellipse-plot
%   SCATTER(X,Y,C,Cov) displays colored ellipses at the locations
%   specified by the vectors X and Y (which must be the same size).  
%   
%   C determines the color-hue of the markers. C should be a vector the
%   same length as X and Y, the values in C are linearly mapped
%   to the colors in the current colormap.
%   
%   Cov determines the ellipse-size and shape each marker. Cov
%   should be a 2 x 2 x nP array with the covariances of the
%   C-values at points [X, Y]. For convenience the input arrays X,
%   Y and C will be converted to column arrays with the (:)
%   operation. In order to avoid sorting confusion it is strongly
%   preferable to arrange these input arrays into column (or row)
%   arrays so that the covariance matrices will be used for the
%   correct points [X, Y, C].
%
%   H = SCATTER_ELLIPSE(...) returns handles to the patch objects created.
% 
% Calling: 
%   h = scatter_ellipse(X,Y,C,Cov[,'covclim',clim,'pclr','r','edgecolor','g'])
% Input:
%   X - double array [nP x 1], x-coordinates of points
%   Y - double array [nP x 1], y-coordinates of points
%   C - double array [nP x 1], value of point, linearly mapped
%       between min and max of current colourmap.
%   COV - covariance matrix for each [Xi, Yi]-point, double array
%         [2 x 2 x nP]
% Optional input arguments (property-value pairs):
%   Covclims - scaling of covariance-ellipse-area (CEA) colour between
%              RGB-triplet (smallest area) towards white (larges
%              area) for color of C. Default value [0.1 1], 1 -
%              RGB-tripplet for the smallest CEA, 0.1 -
%              RGB*0.1+0*[1,1,1] for the point with the largest
%              CEA.
%   pclr - plot-colour, standard matlab-color specification, but
%          with 'rgb' function will plot point with RGB-colour
%          corresponding to C.
%   edgecolor - edgecolor of ellipse, defaults to 'none'.
% Output:
%   h - array with graphics handle to covariance ellipses (plotted
%       with fill)
% Example:
%  X = 12*randn(31,1);
%  Y = 12*randn(31,1);
%  for i1 = 1:31,
%    sx = 1+3*rand(1); 
%    sy = 1+3*rand(1);
%    x = sx*randn(41,1);
%    y = sy*randn(41,1) + randn(1)*x;
%    CvM(:,:,i1) = cov([x,y]);
%    C(i1) = mean(x)*mean(y);
%  end
%  clf
%  h = scatter_ellipse(X,Y,...
%                      medfilt1(C),CvM/3,...
%                      'pclr','rgb',...
%                      'Covclim',[0.15 1],...
%                      'edgecolor','none');
% 
% The idea is that points with larger covariance ellipse should get
% reduced graphical weight while points with smaller covariance
% ellipses should be more clearly visible. The RGB-blending towards
% white could be seen as each point have a fixed amount of pigment
% spread out to color the area of the ellipse.
%  Copyright © Bjorn Gustavsson 20190410, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later
dOPS.Covclim = [0.1 1];
dOPS.pclr  = 'rgb';
dOPS.edgecolor  = 'none';
for i1 = 1:2:numel(varargin)
  curr_option = varargin{i1};
  switch lower(varargin{i1})
   case 'covclim'
    dOPS.Covclim = varargin{i1+1};
   case 'pclr'
    dOPS.pclr = varargin{i1+1};
   case 'edgecolor'
    dOPS.edgecolor = varargin{i1+1};
   otherwise
  end
end
dOPS.Covclim(2) = max(0,min(1,dOPS.Covclim(2)));
dOPS.Covclim(1) = max(0,min(1,dOPS.Covclim(1)));
X = X(:);
Y = Y(:);
C = C(:);
Cmap = colormap;
nCmap = size(Cmap,1);
for i3 = numel(X):-1:1,
  [Xe(i3,:),Ye(i3,:),Aellipse(i3)] = ellipse_xy(Cov(:,:,i3));
end
[Aellipse,iS] = sort(Aellipse);
Xe = Xe(iS,:);
Ye = Ye(iS,:);
X = X(iS);
Y = Y(iS);
C = C(iS);
if max(C) ~= min(C)
  rgbCp = interp1(linspace(0,1,nCmap),...
                 Cmap,...
                 (C-min(C))/(max(C)-min(C)),...
                 'pchip');
else
  rgbCp = repmat(Cmap(end,:),size(C));
end
AeMax = max(Aellipse);
Aemin = min(Aellipse);
Covclim = dOPS.Covclim;
for i1 = numel(X):-1:1,
  resaturation = Covclim(2) - (Aellipse(i1)-Aemin)/(AeMax-Aemin)*diff(Covclim);
  rgbC(i1,:) = rgbCp(i1,:)*resaturation + ...
                   [1 1 1]*(1-resaturation);
end
hold_state = get(gca,'nextplot');
hold on
for i1 = numel(X):-1:1;
  h(i1) = fill(X(i1)+Xe(i1,:),Y(i1)+Ye(i1,:),rgbC(i1,:));
  if strcmp(dOPS.pclr,'rgb')
    plot(X(i1),Y(i1),'.','color',rgbCp(i1,:))
  else
    plot(X(i1),Y(i1),'.','color',dOPS.pclr)
  end
end
set(h,'edgecolor',dOPS.edgecolor)
set(gca,'nextplot',hold_state);


function [x,y,Aellipse] = ellipse_xy(C)
p = linspace(0,2*pi,361); % one-degree steps around 2*pi
[eigvec,lambda] = eig(C); % Eigenvalues and eigen-vectors
xy = [cos(p'),sin(p')] * sqrt(lambda) * eigvec'; % Transformation
x = xy(:,1);
y = xy(:,2);
Aellipse = pi*prod(diag(lambda));
function [xin, yin, y_error] = processVecs(x, y, xrange, y_error)

if size(y, 2) == 1
    x(isnan(y)) = [];
    y(isnan(y)) = [];
end

if nargin > 3
    y_error(isnan(y)) = [];
end

%interpolate the data to be able to fit and get CIs. otherwise the fit is
%impossible since #params ~ #data points.
% scale_interp = 2;
% xq = (min(x): mean(diff(x))/scale_interp : max(x) )';
% yq = interp1(x,y,xq);
xq = x;
yq = y;

%restrict the range of the fit to a monotonically increasing region

if isnan(xrange(1))
    xrange(1) = xq(1);
end
if isnan(xrange(2))
    xrange(2) = xq(end);
end

% x1_ind = find(xq==xrange(1));
% x2_ind = find(xq==xrange(2));
% 
% xin = xq(x1_ind:x2_ind);
% yin =  yq(x1_ind:x2_ind);
xin = xq;
yin = yq;

end
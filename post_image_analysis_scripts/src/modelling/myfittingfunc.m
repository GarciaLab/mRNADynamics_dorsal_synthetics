function forwardModel = myfittingfunc()


 [~, resultsFolder] = getDorsalFolders;
    load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
    
x = dorsalResultsDatabase.dorsalFluoBins( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, '1Dg11') );
y = dorsalResultsDatabase.meanFracFluoEmbryo( ...
       strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, '1Dg11' ) );
y_error = dorsalResultsDatabase.seFracFluoEmbryo( ...
       strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, '1Dg11') );
       
%remove nans from the data for fitting
y_error(isnan(y)) = [];
x(isnan(y)) = [];
y(isnan(y)) = [];

%interpolate the data to be able to fit and get CIs. otherwise the fit is
%impossible since #params ~ #data points. 
scale_interp = 2;
xq = (min(x): mean(diff(x))/scale_interp : max(x) )';
yq = interp1(x,y,xq);

%restrict the range of the fit to a monotonically increasing region
xrange(1,:)  = [1000, 2250];
k = 1;
if isnan(xrange(k, 1))
    xrange(k, 1) = xq(1);
end
if isnan(xrange(k, 2))
    xrange(k, 2) = xq(end);
end
    
x1_ind = find(xq==xrange(k, 1));
x2_ind = find(xq==xrange(k, 2));
    
t = xq(x1_ind:x2_ind);
Y=  yq(x1_ind:x2_ind);


model = @(A, k1, k2, k3, k4, r1, r2, r3, r4, offset) ...
    ... 
  A.* ( (x.*k2.*k3.*(k1+k4+x.*r2)+k3.*(k4+x.*r2).*r3+(k4+r1+x.*r2).*r3.* ...
  r4).*(x.^2.*k2.*(k3+r1).*r2+(k1+r3).*(k3.*k4+(k4+r1).*r4)+x.*(k2.* ...
  k3.*k4+k1.*k2.*(k3+k4+r1)+(k3+r1).*r2.*r3+r2.*(r1+r3).*r4)).^(-1) ) + offset - y;

modelStr = func2str(model);
modelStr = strrep(modelStr, '@(A,k1,k2,k3,k4,r1,r2,r3,r4,offset)', '@(m)');
modelStr = strrep(modelStr, 'x', 't');
modelStr = strrep(modelStr, 'A', 'm(1)');
modelStr = strrep(modelStr, 'k1', 'm(2)');
modelStr = strrep(modelStr, 'k2', 'm(3)');
modelStr = strrep(modelStr, 'k3', 'm(4)');
modelStr = strrep(modelStr, 'k4', 'm(5)');
modelStr = strrep(modelStr, 'r1', 'm(6)');
modelStr = strrep(modelStr, 'r2', 'm(7)');
modelStr = strrep(modelStr, 'k3', 'm(8)');
modelStr = strrep(modelStr, 'k4', 'm(9)');
modelStr = strrep(modelStr, 'offset', 'm(10)');
forwardModel = str2func(modelStr);

end
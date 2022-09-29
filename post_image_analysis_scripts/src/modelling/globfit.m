
% enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3'};

scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73]';

%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = nan(length(enhancers), 2);
xrange(1,:)  = [1000, 2250];  %1Dg
xrange(2, :) = [500, 1750]; %upper limit for %1DgS
xrange(3, 1) = 500; %lower limit for 1DgW
xrange(4, 2) = 1500; %upper limit for %1DgAW


  x1 = dorsalResultsDatabase.dorsalFluoBins( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{1}) );
    y1 = dorsalResultsDatabase.meanFracFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{1} ) );
  x2 = dorsalResultsDatabase.dorsalFluoBins( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{2}) );
    y2 = dorsalResultsDatabase.meanFracFluoEmbryo( ...
           strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{2} ) );
    


% Load input vectors for curves and goal shape for curves
% load inputVec1 inputVec2 goalVec1 goalVec2
% KnownParam = [1, 2, 3];
knownParam = [];

%rate, kd, pkd, offset, omegadp
    y_max = nanmax(y1(:));
    x_max = max(x1);
    %specific for fraction active fitting
    p0 = [y_max; x_max/2 ; 1 ;0; 1];
    lb = [0; 1000; 0; 0; 1];
    ub = [min(1, y_max*2); Inf; Inf; y_max*10; Inf];

[inputVec1 goalVec1] = processVecs(x1, y2, xrange(1, :));
[inputVec2 goalVec2] = processVecs(x2, y2, xrange(2,:));

% How you probably are doing it now
% Setup function handles for individual problems, which return a scalar value that represents the badness of fit
funcFitQual1 = @(unknownParam) calcNormErr( goalVec1, calcCurve( inputVec1, knownParam, unknownParam))
funcFitQual2 = @(unknownParam) calcNormErr( goalVec2, calcCurve( inputVec2, knownParam, unknownParam))
% solve problems individually
fitSol1 = fminsearch(@(x) funcFitQual1 (x), p0)
fitSol2 = fminsearch(@(x) funcFitQual2 (x), p0)

% How you can also do it
% try to solve both problems by adding up their badness of fit
% fitSolBoth= fminsearch(@(x) sum(funcFitQual1(x), funcFitQual2(x)), p0)
fitSolBoth= fminsearch(@(x) funcFitQual1(x)+funcFitQual2(x), p0)


% There are no linear constraints, so set those arguments to [].
A = [];
b = [];
Aeq = [];
beq = [];

fitSolBothCon = fmincon(@(x) funcFitQual1(x)+funcFitQual2(x),p0,A,b,Aeq,beq,lb,ub)


%% Inline functions
function Output = calcCurve(inputVec, knownParam, unknownParam)
    % Some Function which you are trying to fit to data
    model = dorsalFitFunction('simpleWithPol');
    Output = model(unknownParam, inputVec);
end

function errSumAbsNorm = calcNormErr(goalVec, curveVec)
    err = goalVec - curveVec;
    errNorm = err ./ goalVec;
    errSumAbsNorm = sum(abs(errNorm));
end
function [fit, model] = fitDorsalActivity(dlfluobins, activity, varargin)
%possible model types- {'hill', 'simpleWithPol', 'mcwNoPol'}
arguments
    dlfluobins double
    activity double
end
arguments(Repeating)
    varargin
end

modelType = 'hill';
% xScale =10^-(round(log10(max(x(:)))));
% yScale = 10^-(round(log10(max(y(:)))));
xScale = 1;
yScale = 1;
y = activity(:);
x = dlfluobins(:);
x(isnan(y)) = [];
y(isnan(y)) = [];

%scale the problem for better fitting
x = x.*xScale;
y = y.*yScale;
%p(1)=rate coefficient, p(2)=kd, p(3)=hill coefficient p(4) y offset


%rate, kd, hill, y offset
p0 = [max(y) ; max(x)/2 ; 1 ; 0];
lb = [0; 1000.*xScale; 2; -max(y)];
ub = [max(y)*2; Inf; 6; max(y)*10];

%eg for accessing use nlargs{'KD', 'ub'}
nlargs = table(p0, lb, ub, 'RowNames', {'rate', 'KD', 'hill', 'y offset'},...
    'VariableNames', {'p0', 'lb', 'ub'});

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'fixKD')
        nlargs{'KD', : } = varargin{k+1};
    elseif strcmpi(varargin{k}, 'fixRate')
        nlargs{'rate', : } = varargin{k+1};
    elseif strcmpi(varargin{k}, 'fixYOffset')
        nlargs{'y offset', : } = varargin{k+1};
    elseif strcmpi(varargin{k}, 'fixHill')
        nlargs{'hill', : } = varargin{k+1};
    elseif strcmpi(varargin{k}, 'fixW')
        nlargs{'w', :} = varargin{k+1};
    elseif strcmpi(varargin{k}, 'modelType')
        modelType = varargin{k+1};
    end
end

switch modelType
    
    case 'hill'
        %nothing extra to do
    case 'simpleWithPol'
        %add a fifth parameter. this is w = num pol/kd pol
        nlargs{'w', : } = [1, 0, Inf];
        %no y offset for this model
        nlargs{'y offset', :} = [0 0 0]; 
    case 'mwcNoPol'
        nlargs{'y offset', :} = [0 0 0]; 
    otherwise
        error('no valid modeltype')
end

model = dorsalFitFunction(modelType);

options = optimoptions(@lsqcurvefit, 'MaxFunctionEvaluations',...
    size(nlargs, 1)*1000, 'MaxIterations', 4000,...
    'OptimalityTolerance',1E-10,'FunctionTolerance',1E-10,...
    'Display','none');

fit = lsqcurvefit(model, nlargs.p0, x, y, nlargs.lb, nlargs.ub, options);

%gonna transpose here for compatibility with
%other functions
fit = fit';

%rescale
fit(2) = fit(2)./xScale; %kd (x when y is half the plateau)
fit(1) = fit(1)./yScale; %rate coefficient (plateau)
fit(4) = fit(4)./yScale; % y offset

% if length(fit)>4
%     fit(5) = fit(5)./xScale; %kd (x when y is half the plateau)
% end



end
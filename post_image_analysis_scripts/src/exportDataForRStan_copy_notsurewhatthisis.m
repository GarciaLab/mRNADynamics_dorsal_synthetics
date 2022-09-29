function exportDataForRStan(varargin)

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')

expmnt = "affinities"; %affinities, phases or phaff
md = "simpleweak"; %simpleweak, simpleweakdimer, repression, tfdriven, artifact, fourstate
metric = "fluo"; %fraction, fluo
lsq = false;
noOff = true;
nSimu = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
minKD = 200;
maxKD = 1E4;
minw = 1E-2; %1E-2
maxw = 1E1; %1E2
minR = 10;
maxR = 1E3;
displayFigures = true;
wb = true;
fixedKD = NaN; %if this value isn't nan, this value will determine the fixed KD parameter and KD won't be fitted
fixedOffset = NaN;
fixedR = NaN;
fixedw = NaN;
enhancerSubset = {};
scoreSubset = [];
positionSubset = [];
useBatches = false; %fit all the data across embryos and not just means 
noAverage = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if expmnt == "phaff"
    lsq = false;
    noOff = true;
end

if noAverage
    useBatches = false;
end

enhancers_1dg = {'1Dg11'};
enhancers_aff =  {'1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
enhancers_ph = {'1Dg-5', '1Dg-8D'};
scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
positions = [0, -5, -8]';

if strcmpi(expmnt, 'affinities')
    enhancers = [enhancers_1dg, enhancers_aff];
elseif strcmpi(expmnt, 'phases')
    enhancers = [enhancers_1dg, enhancers_ph];
elseif expmnt=="phaff"
    enhancers = [enhancers_1dg, enhancers_aff, enhancers_ph];
end

if ~isempty(enhancerSubset)
    enhancers = enhancerSubset;
    scores= scoreSubset;
    positions = positionSubset;
end


%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = getXRange(enhancers, expmnt);

nSets = length(enhancers);
xo = {};
yo = {};
xs = {};
ys = {};
dsid = [];
T = [];
Y = [];
T_batch = [];
Y_batch = [];
dsid_batch = [];
for k = 1:nSets
    cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k});
    xo{k} = dorsalResultsDatabase.dorsalFluoBins(cond);
    if metric == "fraction"
        yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
        yo_batch{k} = dorsalResultsDatabase.fracFluoEmbryo(cond, :);
    elseif metric == "fluo"
        yo{k} = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond);
        yo_batch{k} = dorsalResultsDatabase.allMaxFluoEmbryo(cond, :);
    end
    
     xo_batch{k} = repmat(xo{k}, [1, size(yo_batch{k}, 2)]);
    [xs_batch{k}, ys_batch{k}]= processVecs(xo_batch{k}, yo_batch{k}, xrange(k, :));
    
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    dsid_batch = [dsid_batch; k*ones(size(xo{k}))];

    if isnan(xrange(k, 1))
        xrange(k, 1) = min(xo{k});
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = max(xo{k});
    end
    
    assert(~isempty(xs{k}));
    
    T = [T; xs{k}];
    Y = [Y; ys{k}];
    T_batch = [T_batch; xo{k}];
    Y_batch = [Y_batch; ys_batch{k}];
end

if ~noAverage
    data.ydata = [T, Y];
    data.dsid = dsid;
    data.X =  [T dsid];
end

data_batch = {};
for k = 1:size(Y_batch, 2)
    data_batch{k}.ydata = [T_batch Y_batch(:, k)];
    data_batch{k}.X =  [T_batch dsid_batch];
    data_batch{k}.dsid = dsid_batch;
end

if useBatches
    disp('Doing batched fits');
    data = data_batch;
end

%%


%rate, kd, hill, y offset
y_max = nanmax(Y(:));
x_max = max(T);


[p0, lb, ub, names] = getInits(expmnt, md, metric ,x_max, y_max, nSets,...
    'minKD', minKD, 'maxKD', maxKD, 'minw', minw, 'maxw', maxw, 'minR', minR, 'maxR', maxR, 'subset', ~isempty(enhancerSubset));

if noOff && metric=="fluo"
    p0 = p0(1:end-1);
    lb = lb(1:end-1);
    ub = ub(1:end-1);
    names = names(1:end-1);
end
enhancers_1dg = {'1Dg11'};
cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers_1dg{1});
% x = dorsalResultsDatabase.dorsalFluoBins(cond);
% Y = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
% x(isnan(Y)) = [];
% Y(isnan(Y)) = [];

x = T;
writetable(table(x, Y), 'C:\Users\owner\Documents\Dorsal-Synthetics-Analysis\dat\dat.csv');
writetable(table(names, p0, lb, ub), 'C:\Users\owner\Documents\Dorsal-Synthetics-Analysis\dat\params.csv');


data2 = {};
datNoAv = table;
    for k = 1:length(enhancers)
        load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'combinedCompiledProjects_allEnhancers')
        ccp = combinedCompiledProjects_allEnhancers(contains({combinedCompiledProjects_allEnhancers.dataSet}, enhancers{k} ) &...
           contains({combinedCompiledProjects_allEnhancers.dataSet}, '2x') &...
            [combinedCompiledProjects_allEnhancers.cycle]==12);
        if metric == "fraction"
            ccp = ccp(cellfun(@any, {ccp.dorsalFluoFeature}));
        end
        isActive = cellfun(@any, {ccp.particleFrames});
        if metric == "fluo"
            
            %make and preprocess the data as a table
            datNoAvTemp = table;
            datNoAvTemp.x = [ccp(isActive).dorsalFluoFeature]';
            datNoAvTemp.Y = cellfun(@max, {ccp(isActive).particleFluo3Slice})';
            datNoAvTemp = rmmissing(datNoAvTemp);
            datNoAvTemp = unique(sortrows(datNoAvTemp));
            datNoAvTemp.dsid = k*ones(size(datNoAvTemp, 1), 1);
            
            
            data2{k}.ydata = [ [ccp(isActive).dorsalFluoFeature]' , cellfun(@max, {ccp(isActive).particleFluo3Slice})'];
        elseif metric == "fraction"
            data2{k}.ydata = [ [ccp.dorsalFluoFeature]' , isActive'];
        end
        simpleWeakOptions.onedsid = true;
        datNoAv = vertcat(datNoAv,datNoAvTemp);
    end
    
writetable(datNoAv, 'C:\Users\owner\Documents\Dorsal-Synthetics-Analysis\dat\datNoAv.csv');    
    1
    
    
end


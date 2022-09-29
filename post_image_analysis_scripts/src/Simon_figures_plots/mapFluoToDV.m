function outStruct = mapFluoToDV(inStruct)

%% load the Reeves data
ReevesDataPath = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox\manuscript\rebuttal_figures\DVaxis_mapping';
Xvalues = load([ReevesDataPath '\Reeves_2012_nc12_X.mat'],'X'); Xvalues = Xvalues.X;
Yvalues = load([ReevesDataPath '\Reeves_2012_nc12_Y.mat'],'Y'); Yvalues = Yvalues.Y;
%plot(Xvalues,Yvalues)

%% fit data to a gaussian and interpolate
f = fit(Xvalues.',Yvalues.','gauss2');
%plot(f,Xvalues,Yvalues)
fittedY = f(linspace(-0.66,0.66,200));
%plot(linspace(-0.66,0.66,200),fittedY,'ko')

fittedY2 = mean([flip(abs(fittedY(1:100))),fittedY(101:200)],2);
%plot(linspace(0,0.66,100),fittedY2,'b*')

%% grab the struct with all data and add DV info
% 
% load('S:\Simon\Dropbox\DorsalSyntheticsDropbox\dorsalResultsDatabase.mat')
% DorsalFluoValues = [combinedCompiledProjects_allEnhancers.dorsalFluoFeature];
% histogram(DorsalFluoValues)

%% rescale Reeves data to the values observed in our experiments
minFluoData = 35;
maxFluoData = 3000;
dynRange = maxFluoData/minFluoData;

maxFluoDV = max(fittedY2);
minFluoDV = min(fittedY2);

ReevesOffset = (maxFluoDV-(dynRange*minFluoDV))/-(dynRange-1);
fittedY3 = (fittedY2-ReevesOffset);%.*dynRange;
fittedY4 = fittedY3*(maxFluoData/max(fittedY3));
%plot(linspace(0,0.66,100),fittedY4,'r.')
AbsPositions = linspace(0,0.66,100);
%% add absolute DV position to the combinceCompiled... struct

for i = 1:length(inStruct)
    i/length(inStruct);
    nucleusFluo = inStruct(i).dorsalFluoFeature;
    if ~isempty(nucleusFluo)
        [dummy,pos] = min(abs(fittedY4-nucleusFluo));
        inStruct(i).AbsDV = AbsPositions(pos);
    else
        inStruct(i).AbsDV = [];
    end
end

outStruct = inStruct;

% %% save results
% clearvars -except combinedCompiledProjects_allEnhancers dataTypes dorsalResults dorsalResultsClean dorsalResultsDatabase field...
% flyregScores names 
% savePath = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox';
% save([savePath '\dorsalResultsDatabase.mat'])




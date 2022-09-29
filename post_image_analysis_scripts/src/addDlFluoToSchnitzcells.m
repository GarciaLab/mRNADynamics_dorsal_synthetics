function schnitzcells = addDlFluoToSchnitzcells(Prefix)

displayFigures = false;

if displayFigures
%     tileFig = figure(); %commented out because it's probably not useful
    holdFig = figure();
end

liveExperiment = LiveExperiment(Prefix);
schnitzcells = getSchnitzcells(liveExperiment);
schnitzcellsOld = schnitzcells;
ncFrames = [zeros(1, 8), liveExperiment.anaphaseFrames'];
FrameInfo = getFrameInfo(liveExperiment);
%FrameNCs = [FrameInfo.nc];
imSize = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine];
resultsFolder = liveExperiment.resultsFolder;

fluoFeatures = []; % SA: what is this for??
% SA: find in which cycle each schnitz mid frame exists to add nc info to
% schnitzcells
for s = 1:length(schnitzcells)
    midFrame = ceil(length(schnitzcells(s).frames)/2); %the middle frame of this schnitz existence
    
    % make sure the schnitz doesn't exist through mitoses
    firstFrame = schnitzcells(s).frames; firstFrame = firstFrame(1);
    lastFrame = schnitzcells(s).frames; lastFrame = lastFrame(end);
    
%     if firstFrame < length(FrameNCs) && length(schnitzcells(s).frames)>1
%         assert(FrameNCs(firstFrame+1)==FrameNCs(lastFrame),'schnitz spanning nuclear cycles found');
%         % if this assert fails fix the call to breakUpSchnitzesAtMitoses
%     end
    
    dif = double(schnitzcells(s).frames(midFrame)) - ncFrames;
    cycle = find(dif>0, 1, 'last' );
    schnitzcells(s).cycle = uint8(cycle);
end

% add the absolute time since the previous mitosis
schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, ncFrames);

%filter out schnitz that last too short or that are too close to the edge
schnitzcells = filterSchnitz(schnitzcells, imSize);


%filter out schnitz that have an empty fluo field.
   schnitzcells = filterSchnitzFurther(...
    schnitzcells);

approvedSchnitzes = find([schnitzcells.Approved]);

for s = approvedSchnitzes
    
    concavityThreshold = 0.5; %this is the quadratic coefficient of a parabola that we fit to
    %to the curve of fluo over to see if it's is pointing up or down.
    schnitzcells(s).FluoTimeTrace = single(ExtractDlFluo(schnitzcells(s).Fluo,concavityThreshold));
    
    if schnitzcells(s).cycle ~= 14
        midCycle = floor((ncFrames(schnitzcells(s).cycle) + ncFrames(schnitzcells(s).cycle+1))/2);
    else
        midCycle = 1;
    end
    midCycleFrame = find(schnitzcells(s).frames==midCycle);
    schnitzcells(s).midCycleFrame = find(schnitzcells(s).frames==midCycle);
    
    schnitzcells(s).FluoTimeTraceSmooth = smooth(single(ExtractDlFluo(schnitzcells(s).Fluo,concavityThreshold)), 5);
    
%     % SA: this doesn't get used, I don't see the purpose
%     if schnitzcells(s).cycle ~= 14
%         midCycleSmooth = floor((ncFrames(schnitzcells(s).cycle) + ncFrames(schnitzcells(s).cycle+1))/2);
%     else
%         midCycleSmooth = 1;
%     end
%     midCycleFrameSmooth = find(schnitzcells(s).frames==midCycleSmooth);
    
    
    %different characterizations of intensity
    schnitzcells(s).fluoMid= single(schnitzcells(s).FluoTimeTrace(midCycleFrame)); %middle fluorescence
    schnitzcells(s).fluoMidSmooth = single(schnitzcells(s).FluoTimeTraceSmooth(midCycleFrame)); %middle fluorescence smoothed
    
    schnitzcells(s).FluoFeature =  schnitzcells(s).fluoMid; %the preferred intensity
    
    % SA: why would this happen? maybe when the schnitz doesn't exist in
    % the mid cycle frame? in that case we can take the fluo from the
    % neighbors
    if isempty(schnitzcells(s).FluoFeature)
        schnitzcells(s).FluoFeature = NaN;
    end
    
    schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);
    
    
    fluoFeatures = [fluoFeatures, schnitzcells(s).FluoFeature]; %#ok<AGROW>
    
    if displayFigures && schnitzcells(s).cycle == 12 && ~isnan(schnitzcells(s).FluoFeature)
        
%         figure(tileFig)
%         nexttile;
%         plot(schnitzcells(s).frames, schnitzcells(s).FluoTimeTrace, '-k');
%         hold on
%         plot(midCycle, schnitzcells(s).FluoFeature, 'ob');
%         hold off
%         xticks([]);
%         yticks([round(min(schnitzcells(s).FluoTimeTrace)), round(max(schnitzcells(s).FluoTimeTrace))]);
        %                 ylim([0, 3000]);
        
        figure(holdFig)
%         plot(1:length(schnitzcells(s).FluoTimeTrace), schnitzcells(s).FluoTimeTrace, '-k');
        plot(schnitzcells(s).timeSinceAnaphase, schnitzcells(s).FluoTimeTrace, '-k');
        hold on
%         plot(midCycleFrame, schnitzcells(s).FluoFeature, 'ob');
        ylim([0, 3000]);
        
    end
    
    
end

% SA: make sure we have a reasonable number of nc12 nuclei with fluo
nc12Schnitz = [schnitzcells.cycle] == 12;
fluoFeatureNc12Schnitz = [schnitzcells(nc12Schnitz).FluoFeature];
usefulNc12Schnitz = sum(~isnan(fluoFeatureNc12Schnitz));
assert(usefulNc12Schnitz >= 11,['very few good nuclei in ' Prefix])




save([resultsFolder,filesep,Prefix,'_lin.mat'], 'schnitzcells');


end
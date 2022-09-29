function schnitzcells =...
    integrateSchnitzFluo_clusters(Prefix, schnitzcells, FrameInfo)

saveFlag = false;

liveExperiment = LiveExperiment(Prefix);

Ellipses = getEllipses(LiveExperiment(Prefix))

if nargin == 1
        
    DropboxFolder = liveExperiment.userResultsFolder;
    Channels = liveExperiment.Channels;
    
    FrameInfo = getFrameInfo(liveExperiment); 
    
    schnitzcells = getSchnitzcells(liveExperiment);

    saveFlag = true;
    
end

schnitzPath = [liveExperiment.resultsFolder, filesep, Prefix, '_lin.mat'];

Channels = liveExperiment.Channels;

InputChannelIndexes = find(contains(Channels, 'input', 'IgnoreCase', true));

if isempty(InputChannelIndexes)
    warning(['No input channel found. Check correct definition in MovieDatabase.',...
        ' Input channels should use the :input notation.'])
    return;
end

movieMat = getMovieMat(liveExperiment);

numFrames = length(FrameInfo);


%Create the circle that we'll use as the mask
IntegrationRadius=2;       %Radius of the integration region in um
IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels

if ~mod(IntegrationRadius,2)
    IntegrationRadius=IntegrationRadius+1;
end

Circle=false(3*IntegrationRadius,3*IntegrationRadius);
Circle=double(MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
    1.5*IntegrationRadius+0.5,1));

%Initialize fields
schnitzcells(1).Fluo = [];

if sum(InputChannelIndexes)
    
    %Extract the fluroescence of each schnitz, for each channel,
    %at each time point
    
    %Get the image dimensions
    PixelsPerLine=FrameInfo(1).PixelsPerLine;
    LinesPerFrame=FrameInfo(1).LinesPerFrame;
    %Number of z-slices
    nPadding = 2;
    
    nSlices=FrameInfo(1).NumberSlices + nPadding;
    
    
    %Generate reference frame for edge detection
    refFrame = ones(LinesPerFrame,PixelsPerLine, nSlices);
    convRef = convn(refFrame, Circle, 'same');
    edgeMask = convRef~=sum(Circle(:));
    chIndex = 0;
    
    for ChN=InputChannelIndexes
        chIndex = chIndex+1;
        h=waitbar(0,['Extracting nuclear fluorescence for channel ',num2str(ChN)]);
        try waitbar(CurrentFrame/numFrames,h); catch; end
        
        % SA: I don't like the assumption that nuclear fluorescence is
        % homogeneous. Instead, let's look at the distribution of
        % fluorescence values in nucleus, each slice and each frame. Then
        % we can apply some criterion to leave aggregates/clusters out of
        % the fluorescence calculation.
        % this horrible nested for loop is to create an array called
        % 'movieFluoStack' which has the same dimensions as the raw data
        % but in the position of each ellipse it contains fluorescence
        % values.

        tempSchnitz = schnitzcells;

        for CurrentFrame= 33% 1:numFrames  
            
            getMovieFrame(liveExperiment,CurrentFrame,ChN)
            movieFluoStack = nan(FrameInfo(1).PixelsPerLine,FrameInfo(1).LinesPerFrame,2+FrameInfo(1).NumberSlices,length(FrameInfo));     
        
            radiusScale = 0.75; %in relation to the radius in Ellipses.mat
            ellipseMask = makeNuclearMask(Ellipses{CurrentFrame},[PixelsPerLine LinesPerFrame],radiusScale);
            labeledMask = bwlabel(ellipseMask);
            zStackRaw = getMovieFrame(liveExperiment,CurrentFrame,ChN); %load the raw data
            zStackMask = repmat(labeledMask,1,1,size(zStackRaw,3));
            for e = 1:max(labeledMask(:)) %loop over ellipses in this frame
                singleEllipseMask = labeledMask == e;
                singleEllipseStack = immultiply(repmat(singleEllipseMask,1,1,size(zStackRaw,3)),zStackRaw);
                ellipseSliceFluo = zeros(1,size(zStackRaw,3));
                for z = 2:size(zStackRaw,3)-1 %loop over z slices, first and last are blank
                    ellipseSlice = singleEllipseStack(:,:,z);
                    ellipseSlice(ellipseSlice==0) = nan;
                    %histogram(ellipseSlice)
                    ellipseSliceFluo(z) = nanmean(ellipseSlice);
                    sliceStack = zStackMask(:,:,z); sliceStack(sliceStack==e) = ellipseSliceFluo(z);
                    zStackFluos(:,:,z) = sliceStack;
                end %finish looping over z slices for this ellipse
            end %finish looping over ellipses in this frame
            movieFluoStack(:,:,:,CurrentFrame) = zStackFluos;            
            
            if ~isempty(movieMat)
                imStack = double(movieMat(:,:,:, CurrentFrame, ChN));
            else
                imStack = double(getMovieFrame(liveExperiment, CurrentFrame, ChN));
            end
            % 
            convImage = imfilter(imStack, Circle, 'same');
            convImage(edgeMask) = NaN;
            zStackFluos(edgeMask) = NaN;
            
            for j=1:length(tempSchnitz)
                
                cenx = min(max(1,round(tempSchnitz(j).cenx(tempSchnitz(j).frames==CurrentFrame))),PixelsPerLine);
                ceny = min(max(1,round(tempSchnitz(j).ceny(tempSchnitz(j).frames==CurrentFrame))),LinesPerFrame);
                tempSchnitz(j).Fluo_noClusters(tempSchnitz(j).frames==CurrentFrame,:,chIndex) = single(zStackFluos(ceny,cenx,:));
                tempSchnitz(j).Fluo(tempSchnitz(j).frames==CurrentFrame,:,chIndex) = single(convImage(ceny,cenx,:));
    
                    
            end %loop of nuclei in a frame
            
        end %loop over frames
        
        schnitzcells = tempSchnitz;
        
        
        try close(h); end
    
    end %loop over channels
    
else 
    
    error('Input channel not recognized. Check correct definition in MovieDatabase.Input channels should use the :input notation.');

end 

if saveFlag
    if whos(var2str(schnitzcells)).bytes < 2E9
        save(schnitzPath, 'schnitzcells', '-v6');
    else
        save(schnitzPath, 'schnitzcells', '-v7.3', '-nocompression');
    end
    
end

end
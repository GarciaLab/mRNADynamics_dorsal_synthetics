function Spots =...
    ...
    addSpot2(prefix)

liveExperiment = LiveExperiment(prefix);
FrameInfo = getFrameInfo(liveExperiment);

segmentChannel = 1;
CurrentChannel = 2;

use_integral_center = 1;


Spots_in = getSpots(liveExperiment);

spots1 = Spots_in{segmentChannel};

numFrames = length(spots1);

Spots = repmat(struct('Fits', []), 1, numFrames);
Spots = {spots1, Spots};


for CurrentFrame = 1:numFrames
    
    %report progress every fifth frame
    if ~mod(CurrentFrame, 5), disp(['Segmenting frame ',...
            num2str(CurrentFrame), '...']); end
    
    spots1Frame = spots1(CurrentFrame).Fits;
    
    imStack = getMovieFrame(liveExperiment, CurrentFrame, CurrentChannel);
    
    for spotIndex = 1:length(spots1Frame)
        
        spot1 = spots1(CurrentFrame).Fits(spotIndex);
        
        SpotsIndex = length(Spots{CurrentChannel}(CurrentFrame).Fits)+1;
         
        FitCell = cell(1, length(spot1.z));
        
        zIndex = 0;
        
        for z = spot1.z
            
            zIndex = zIndex + 1;
            
            spot1x = spot1.xDoG(zIndex);
            spot1y = spot1.yDoG(zIndex); 
         
        
            ConnectPositionx = double(round(spot1x));
            ConnectPositiony = double(round(spot1y));
            
            spotsIm = imStack(:, :, z);
            try
                imAbove= imStack(:, :, z+1);
                imBelow= imStack(:, :, z-1);
            catch
                imAbove = nan(size(spotsIm,1),size(spotsIm,2));
                imBelow = nan(size(spotsIm,1),size(spotsIm,2));
            end
            
            Threshold = min(min(spotsIm));
            dog = spotsIm;
            dog_thresh = dog >= Threshold;
            [dog_label, ~] = bwlabel(dog_thresh);
            microscope = FrameInfo(1).FileMode;
            show_status = 0;
            fig = [];
            k = 1; %This is supposed to be the index for the particles in an image.
            %However, this image only contains one particle
            neighborhood_px = round(1300 / liveExperiment.pixelSize_nm); %nm
            %Get the information about the spot on this z-slice
            
            Fit = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, dog_label, dog, neighborhood_px, liveExperiment.snippetSize_px, ...
                liveExperiment.pixelSize_nm, show_status, fig, microscope,...
                [1, ConnectPositionx, ConnectPositiony], [ConnectPositionx, ConnectPositiony], '', CurrentFrame, [], z);
            
            if isempty(Fit)
                keyboard
            end
            
            FitCell{zIndex} = Fit;
            
            
        end
        Fits = [];
        
        for zIndex = 1:length(spot1.z)
            if ~isempty(FitCell{zIndex})
                fieldnames = fields(FitCell{zIndex});
                if isempty(Fits)
                    Fits = FitCell{zIndex};
                else
                    for field = 1:length(fieldnames)
                        Fits.(fieldnames{field}) =  [Fits.(fieldnames{field}), FitCell{zIndex}.(fieldnames{field})];
                    end
                end
            end
        end
        %
        % if isempty(Fits)
        %     breakflag = true;
        % end
        
        force_z = spot1.brightestZ;
        
        [Fits,~, ~] = findBrightestZ(Fits, -1, use_integral_center, force_z, []);
        
        
        if SpotsIndex ~= 1
            if ~isempty(setdiff(fields(Spots{CurrentChannel}(CurrentFrame).Fits), fields(Fits)))...
                    | ~isempty(setdiff(fields(Fits),fields(Spots{CurrentChannel}(CurrentFrame).Fits)))
                [Spots{CurrentChannel}(CurrentFrame).Fits, Fits] = addFields(Spots{CurrentChannel}(CurrentFrame).Fits, Fits);
            end
            Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex) = Fits;
        else
            Spots{CurrentChannel}(CurrentFrame).Fits = Fits;
        end
        
    end
end

save([liveExperiment.resultsFolder,...
        filesep, 'Spots.mat'], 'Spots', '-v6');

TrackmRNADynamics(prefix);

end

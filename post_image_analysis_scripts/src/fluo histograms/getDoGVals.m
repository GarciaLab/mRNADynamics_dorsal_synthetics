function vals = getDoGVals(Prefix, shouldMaskNuclei)

if nargin==1
    shouldMaskNuclei = false;
end

liveExperiment = LiveExperiment(Prefix);

DogOutputFolder = [liveExperiment.procFolder, 'dogs', filesep];
dogDir = dir([DogOutputFolder, '*_ch0','*.*']);
dogDir = {dogDir.name};


FrameInfo = getFrameInfo(liveExperiment);
Ellipses = getEllipses(liveExperiment);


nFrames = length(FrameInfo);



vals = [];

for f = 1:nFrames
    
    dog = double(imreadStack([DogOutputFolder, filesep, dogDir{f}]));
    dog = (dog/100) - 100;
%     dog = dog - quantile(dog(:), .1); 
%     dog = dog - min(dog(:)); 
%     dog = a ./ max(a(:));

    if shouldMaskNuclei
        
        nuclearMask = ones(size(dog, 1), size(dog, 2));
        ellipsesFrame = Ellipses{f};
        nuclearMask = makeNuclearMask(ellipsesFrame, [size(dog,1), size(dog,2)], 1);    
        
        %stack the nuclear mask to create cylinders out of the nuclei
        nuclearMask = repmat(nuclearMask, [1 1 size(dog, 3)]);
        
        dog = dog.*nuclearMask;
        
    end
    
    %filter out the z slices away from the center
    for z = 1:size(dog, 3)
        if z<6 || z>12 %take just middle slices
            dog(:, :, z) = 0;
        end
    end
   
%     vals(f) = log10(max(dog(:))+1); %#ok<AGROW>
    vals = [vals, log10(dog(dog > quantile(dog(:), .9))+1)']; %#ok<AGROW>
    %imshow(dog, []);
    
end

% vals(vals==0) = NaN;

end
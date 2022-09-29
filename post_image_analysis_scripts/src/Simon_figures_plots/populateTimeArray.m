function [particleNCFluo] = populateTimeArray(timeArray,particleTime,particleFluo,frameRate)

% timeArray contains the absolute time,
% particleTime contains the time at which each particle fluo observation was
% made counting since anaphase
% particleFluo contains the fluo observations

% particleNCFluo is vector the same size as timeArray
% we want to populate particleNCFluo with either 0 or a fluo value.
% to do this we'll loop over particleTime and ask what's the closest
% absolute time in timeArray.

%% interpolate particle fluorescence
particleNCFluo = zeros(size(timeArray));

if length(particleTime) > 1
    interpPoints = length(timeArray);
    interpTime = linspace(min(particleTime),max(particleTime),interpPoints);
    interpFluo = interp1(particleTime,particleFluo,interpTime);
else
    interpTime = particleTime;
    interpFluo = particleFluo;
    %plot(interpTime/60,interpFluo)
end

for i = 1:length(timeArray)    
    [diff,interpTimeIdx] = min(abs(interpTime-timeArray(i)));   
    if diff < frameRate        
        particleNCFluo(i) = interpFluo(interpTimeIdx);
    end
end


        
    
end

% 
% for i = 1:length(particleTime)
%     
%     time = particleTime(i);
%     [distanceToIntrpTime,Idx1] = min(abs(interpTime-time));
%     [distanceToAbsTime,Idx2] = min(abs(timeArray-time));
%     particleFluoVal = interpFluo(Idx1);
%        
%     if distanceToAbsTime <= (frameRate+1)*3
%         particleNCFluo(Idx2) = particleFluoVal(i);
%     end
% end

    


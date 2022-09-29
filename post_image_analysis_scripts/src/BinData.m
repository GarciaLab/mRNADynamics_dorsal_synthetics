function [binnedvector] = BinData(vector,binValues)
% first argument = vector array to be binned, 
% second arg = vector array containing the upper limit of each of the bins
%%

%Range = Max-Min;
%BinWidth = Range/Bins;
binnedvector =[]; %vector to fill with binning info

for value = 1:length(vector)
    dataValue = vector(value);
    %find the nearest bin limit
    distancesToBinLimits = abs(binValues - dataValue);
    [~,binnedvector(value)] = min(distancesToBinLimits);
    
    
%     if distance == 0; %if the value = the minimum value, then it's in bin 1
%         valuebin = 1;
%     elseif dataValue == NaN
%         valuebin = NaN;
%     else
%         valuebin = ceil(distance/BinWidth);
%     end
    %binnedvector(value) = binValues(index); %fill the output vector
end


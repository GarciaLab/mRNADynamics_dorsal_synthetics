function dorsalResults = addKeysToDorsalResults(dorsalResults)

 dataType = dorsalResults{1}.DataType;
 
 %nc14 is problematic and mostly uninteresting. let's just do 12
    for ncIndex = 1
        
        %this field is problematic. 
         dorsalResults{ncIndex} = rmfield(dorsalResults{ncIndex}, 'allmrnasnomean');

%         try
            setLength = length(dorsalResults{ncIndex}.dorsalFluoBins); 
            dorsalResults{ncIndex} = rmfield(dorsalResults{ncIndex}, 'DataType');
            [dorsalResults{ncIndex}.DataType{1:setLength}] = deal(dataType);
            dorsalResults{ncIndex}.DataType = dorsalResults{ncIndex}.DataType';
            dorsalResults{ncIndex}.nc = ones(setLength, 1)*(ncIndex + 11); %nc offset
            
            %add whether the dataset is 1x or 2x
            if contains(dataType, 'FFF')
               dorsalResults{ncIndex}.mother = cellstr(repmat("FFF", setLength, 1));
            elseif contains(dataType, '2x')
                dorsalResults{ncIndex}.mother = cellstr(repmat("2x", setLength, 1));
            end
  %%          
             %add a dad enhancer attribute
             dataTypeTokens = split(dataType, '_');
             enhancer = dataTypeTokens{1};
             %make a special exception since this name was mistakenly not
             %standardized
             if strcmpi(dataType, '1DG_2xDl')
                 enhancer = [enhancer, '11'];
             end
             dorsalResults{ncIndex}.enhancer = cellstr(repmat(enhancer, setLength, 1));
 %%               
            
            %transpose to conform with the rest of the data
            if size(dorsalResults{ncIndex}.dorsalFluoBins, 1) == 1
                dorsalResults{ncIndex}.dorsalFluoBins = dorsalResults{ncIndex}.dorsalFluoBins';
            end
            
%         catch
%             warning(['skipping: ', num2str(ncIndex)])
%         end
        
    end
    
end
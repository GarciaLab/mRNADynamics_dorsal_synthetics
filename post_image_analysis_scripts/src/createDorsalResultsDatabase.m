function dorsalResultsDatabase =...
    createDorsalResultsDatabase(dataTypes)


if nargin<1
    dataTypes = unique({'1Dg11_2xDl_FFF','1Dg-20(11)_2xDl','1Dg-8D_FFF', '1Dg11_2xDl', '1DgW_2x_Leica',...
        '1DgW_FFF', '1DgVW_FFF', '1Dg11_FFF', '1Dg-5_2xDl','1Dg-5_FFF'...
        '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl', '1DgSVW2_2xDl',...
        '1DgVW_2xDl','TwiPEv5(7)_2xDl','2Dgc_FFF','2Dgc_2xDl','1Dg-12_6_2xDl','1Dg11_2xDl_FFF',...
        '1DgW_2xDl_FFF','1DG_VW_2xDl_FFF'});
end




dorsalResultsDatabase = struct;
combinedCompiledProjects_allEnhancers = [];

for i = 1:length(dataTypes)
    i/length(dataTypes)
    
    [~, resultsFolder] = getDorsalFolders;
    
    dorsalResults = createDorsalResults(dataTypes{i}); 

    dorsalResultsClean = addKeysToDorsalResults(dorsalResults);
    
    %skipping nc14
%     dorsalResultsClean = dorsalResultsClean(1:2);
%and 13 now
    dorsalResultsClean = dorsalResultsClean(1);
    
    for ncIndex = 1:length(dorsalResultsClean)
        
        fields = fieldnames(dorsalResultsClean{ncIndex});
        
        if isempty(fieldnames((dorsalResultsDatabase)))
            emptyCell = cell(length(fields),1);
            dorsalResultsDatabase= cell2struct(emptyCell,fields);
        end
        
        for f = 1:length(fields)
            
            
            if isvector(dorsalResultsClean{ncIndex}.(fields{f}))
                
                dorsalResultsDatabase.(fields{f}) =...
                    [dorsalResultsDatabase.(fields{f})
                    dorsalResultsClean{ncIndex}.(fields{f})];
                
            elseif isnumeric(dorsalResultsClean{ncIndex}.(fields{f}))
                
                %matrices require special care. concacatenate
                %and pad with NaNs since data sets
                %have different numbers of embryos
                dorsalResultsDatabase.(fields{f}) =...
                    padCatMatrices(dorsalResultsDatabase.(fields{f}),...
                    dorsalResultsClean{ncIndex}.(fields{f}) );
                
            end
            
        end
        
    end %nc loop
    try
        %%
        %create a database of all traces for every experiment. 
        load([resultsFolder,filesep,dataTypes{i},filesep,'combinedCompiledProjects.mat'], 'combinedCompiledProjects');
        if i == 1
            combinedCompiledProjects_allEnhancers = combinedCompiledProjects; 
        else
            combinedCompiledProjects_allEnhancers =...
                [combinedCompiledProjects_allEnhancers, combinedCompiledProjects]; %#ok<AGROW>
        end
        clear combinedCompiledProjects; %to ensure we don't accidentally add the same set twice
        %%
    end

end %dataType loop


%%
%now let's add affinity information
names = {'1Dg11', '1DgS2', '1DgAW3', '1DgW', '1DgSVW2', '1DgVW', '1DgVVW3'}';
flyregScores = [6.23, 5.81, 5.13, 5.39, 4.80, 4.29, 4.73]';

dorsalResultsDatabase.patserScore = nan(length(dorsalResultsDatabase.enhancer), 1);
for k = 1:length(names)
    for j = 1:length(dorsalResultsDatabase.enhancer)
        if strcmpi(dorsalResultsDatabase.enhancer{j}, names{k})
            dorsalResultsDatabase.patserScore(j) = flyregScores(k);
        end
    end
end
%%


save([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], '-v6');
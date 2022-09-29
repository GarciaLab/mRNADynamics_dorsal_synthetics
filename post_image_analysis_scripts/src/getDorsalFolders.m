function [dataFolder, resultsFolder, PreProcessedFolder, ProcessedFolder] = getDorsalFolders()

    configValues = csv2cell('ComputerFolders.csv', 'fromfile');
    try
        dataFolder = getConfigValue(configValues, 'DorsalSynthetics');
    catch
        dataFolder = '';
    end
    resultsFolder = getConfigValue(configValues, 'DorsalSyntheticsDropbox');
    PreProcessedFolder = [dataFolder, filesep, 'PreProcessedData'];
    ProcessedFolder = [dataFolder, filesep, 'ProcessedData'];
    
end

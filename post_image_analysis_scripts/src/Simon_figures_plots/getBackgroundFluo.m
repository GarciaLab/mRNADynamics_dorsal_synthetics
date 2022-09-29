function offsetFluo = getBackgroundFluo(PrefixList)

PrefixList = {'2020-10-23-dl1-MCPmCh_VenusBackground_1','2020-10-23-dl1-MCPmCh_VenusBackground_2'};

DynamicsResultsPath = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox';
DataType = 'dl1_MCP_background';
%addDVStuffToSchnitzCells(DataType);

Offsets = [];
for p= 1:length(PrefixList)
    Prefix = PrefixList{p};
    addDlFluoToSchnitzcells(Prefix);
    prefixResultsPath = [DynamicsResultsPath '\' Prefix '\' Prefix '_lin.mat'];
    schnitzcells = load(prefixResultsPath);
    schnitzcells = schnitzcells.schnitzcells;
    schnitzNCs = [schnitzcells.cycle];
    nc12IDx = [1:length(schnitzcells)] .* double(schnitzNCs == 12);
    nc12IDx = nc12IDx(nc12IDx>0);
    
    AllFluoFeatures = [];
    for n = nc12IDx
        schnitzFluoFeature = schnitzcells(n).FluoFeature;
        AllFluoFeatures = [AllFluoFeatures schnitzFluoFeature];
    end 
    Offsets = [Offsets nanmean(AllFluoFeatures)];
end

offsetFluo = nanmean(Offsets);
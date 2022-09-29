% number of nuclei per Dorsal bin normalized to number of embryos that
% include that bin

%% Load stuff
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])%, 'dorsalResultsDatabase')
b = load('S:\Simon\Dropbox\DorsalSyntheticsDropbox\dorsalResultsDatabase.mat');

Nbins = 19; %number of Dorsal concentration bins

% make a substruct containing only nc12 nuclei
AllNc12Struct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);

NucleiPerFluoBin = zeros(1,19);
EmbryosPerFluoBin = zeros(1,19);

%% Loop over fluo bins to count number of nuclei per bin

for nuc = 1:length(AllNc12Struct)
    nucleusFluoBin = AllNc12Struct(nuc).dorsalFluoBin;
    if ~isempty(nucleusFluoBin) && ~isnan(nucleusFluoBin)
        NucleiPerFluoBin(nucleusFluoBin) = NucleiPerFluoBin(nucleusFluoBin)+1;
    end
end
plot(linspace(0,4500,19),NucleiPerFluoBin,'-ko','LineWidth',2)
ylabel('number of nuclei')
xlabel('Dorsal fluorescence (AU)')
%% Loop over embryos (prefixes) to
Prefixes = unique({AllNc12Struct.prefix});

for p = 1:length(Prefixes)
    %p/length(Prefixes)
    PrefixName = Prefixes{p};
    % make substruct containing just this prefix
    EnhancerStruct = AllNc12Struct(strcmpi({AllNc12Struct.prefix}, PrefixName));
    thisEmbryoBins = [EnhancerStruct.dorsalFluoBin];
    thisEmbryoBins = thisEmbryoBins(~isempty(thisEmbryoBins));
    thisEmbryoBins = thisEmbryoBins(~isnan(thisEmbryoBins));
    EmbryosPerFluoBin(thisEmbryoBins) = EmbryosPerFluoBin(thisEmbryoBins)+1;
end

plot(linspace(0,4500,19),EmbryosPerFluoBin,'-ro','LineWidth',2)
ylabel('number of embryos that cover this bin')
xlabel('Dorsal fluorescence (AU)')

%%
plot(linspace(0,4500,19),NucleiPerFluoBin./EmbryosPerFluoBin,'-bo','LineWidth',2)
ylabel('number of nuclei per bin / number of embryos that cover bin')
xlabel('Dorsal fluorescence (AU)')
ylim([15 45])
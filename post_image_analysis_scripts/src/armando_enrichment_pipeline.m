
datasets = {'1Dg-12_6_2xDl', '1Dg-20(11)_2xDl','1Dg-8D_FFF', '1DgW_2x_Leica',...
    '1DgW_FFF', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF', '1Dg-5_2xDl',...
    '1DgS2_2xDl', '1DgAW3_2xDl', '1DgVVW3_2xDl', '1Dg-8D_2xDl', '1DgSVW2_2xDl',...
    '1DgVW_2xDl', '1Dg11_2xDl','TwiPEv5(7)_2xDl'};

names0 = {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'}';
scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';


scores = [6.23, 5.81, 5.13, 4.73, 6.23, 4.80, 4.29, 6.23, 8];
positions = [-5, 0, 0, 0, -8, 0, 0, 0, 0];

% dataType = '1Dg11_2xDl';
% dataType = '1DgS2_2xDl';
dataType = '1DgVW_2xDl';
% dataType = 'TwiPEv5(7)_2xDl'

rmaxes = [];
rstes = [];
names = [];
for k = 1:length(datasets)
  
    dataType = datasets{k};
   
% startup

% addpath(genpath('S:\Armando\nick_enrichment_fork\tf_enrichment_pipeline'));
try
% main01_compile_traces(dataType,'firstNC', 12)

% main02_sample_local_protein(dataType,'ignoreQC', true, 'max_nucleus_radius_um', 6,'segmentNuclei', 1, 'NumWorkers', 24, 'use3DSpotInfo', false);
% main02_sample_local_protein(dataType,'ignoreQC', true, 'max_nucleus_radius_um', 6,'segmentNuclei', 1, 'NumWorkers', 1, 'displayFigures', true);
% main02_sample_local_protein(dataType,'ignoreQC', true, 'max_nucleus_radius_um', 6,'segmentNuclei', 1, 'NumWorkers', 24, 'segmentationMethod', 2);

main04_make_exploratory_figs(dataType)
load(['S:\Armando\Dropbox\LocalEnrichmentFigures\PipelineOutput\', dataType, '\basic_figs\radial_enrichment_data.mat']);
rmaxes = [rmaxes, radial_enrichment_data.r_spot_mean(1) / radial_enrichment_data.r_control_mean(1)];
rstes = [rstes, radial_enrichment_data.r_spot_ste(1) / radial_enrichment_data.r_control_mean(1)];
names = [names, string(dataType)];

catch
    continue
end

end
dat = table(names',rmaxes', rstes', scores', positions');
dat = sortrows(dat, 4);
toDelete = dat.Var1 == "1Dg-5_2xDl" | dat.Var1 == "1Dg-8D_2xDl" | dat.Var1 == "TwiPEv5(7)_2xDl";
dat(toDelete, :) = [];
figure;
errorbar(dat.Var4, dat.Var2, dat.Var3, '.');
ylabel('enrichment')
xlabel('affinity (patser score)')
ylim([1.03, 1.06]);

figure;
dat = sortrows(dat, 5);
errorbar(dat.Var5, dat.Var2, dat.Var3, '.');
ylabel('enrichment')
xlabel('position relative to 1Dg (bp)')


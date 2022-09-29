function plotAllTraces()

[~, resultsFolder] = getDorsalFolders;

load(resultsFolder+"dorsalResultsDatabase.mat", 'combinedCompiledProjects_allEnhancers')

dat = combinedCompiledProjects_allEnhancers; 

%perform some initial QC 
dat = dat([dat.cycle] == 12);
dat = dat(~isnan([dat.dorsalFluoFeature]));
dat = dat(~cellfun(@isempty, {dat.particleFluo3Slice}));
dat = dat(~cellfun(@(x, y) contains(x, '1Dg11'), {dat.dataSet}));
dat = dat(cellfun(@length, {dat.particleFluo3Slice}) > 10); %10 is pretty arbitrary. just testing 
%%
figure;
for k = 1:length(dat)
    
    if rand > .97
        plot(dat(k).particleTimeSinceAnaphase,dat(k).particleFluo3Slice, 'LineWidth', 1)
        hold on
    end

end

title('Random 5% of particle fluorescence traces from DBS_6.29')
ylabel('Particle fluorescence (a.u.)')
xlabel('Time into nc12 (min)')
ylim([0, 600])
xlim([3, 8])


end
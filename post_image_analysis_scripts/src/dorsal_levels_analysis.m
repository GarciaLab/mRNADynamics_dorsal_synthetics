Prefix = '2020-09-17-en2DlV_en2mcpMch_hisirfp_HET_1'; %very ventral movie
schnitzcells = getSchnitzcells(LiveExperiment(Prefix));

figure(1);
ch = 1;
zdim = 2;
for k = 1:length(schnitzcells)
    plot(schnitzcells(k).frames, nanmax(schnitzcells(k).Fluo(:, :, ch),[],zdim))
    hold on
end

figure(2);
ch = 2;
for k = 1:length(schnitzcells)
    plot(schnitzcells(k).frames, nanmax(schnitzcells(k).Fluo(:, :, ch),[],zdim))
    hold on
end

Prefix = '2020-07-01-1DgAW3_2xDl_3_mcpQuant'
schnitzcells = getSchnitzcells(LiveExperiment(Prefix));

figure(3);
ch = 1;
zdim = 2;
for k = 1:length(schnitzcells)
    plot(schnitzcells(k).frames, nanmax(schnitzcells(k).Fluo(:, :, ch),[],zdim))
    hold on
end

figure(4);
ch = 2;
for k = 1:length(schnitzcells)
    plot(schnitzcells(k).frames, nanmax(schnitzcells(k).Fluo(:, :, ch),[],zdim))
    hold on
end
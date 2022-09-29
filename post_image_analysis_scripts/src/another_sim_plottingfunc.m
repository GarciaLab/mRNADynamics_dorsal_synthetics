%%
load('C:\Users\owner\Dropbox\DorsalSyntheticsDropbox\simulations\tf_paramsearch_entry_exitOnlyOnOffStates_.mat')
% load('C:\Users\owner\Dropbox\DorsalSyntheticsDropbox\simulations\archive\tf_paramsearch_entryexit_exitOnlyOnOffStates_.mat')

if isfield(params, 'pi1s')
    params.pi_exits = params.pi1s;
end
if isfield(params, 'pi2s')
    params.pi_entries = params.pi2s;
end

nPlots = length(params.cs);
ni = @(v, x) nearestIndex(v, x);

kd = ni(params.kds, 1E3);
kd2 = ni(params.kds, 1E4);
kd3 = ni(params.kds, 1E5);

if params.model == "entry"
    %this is a "good" solution with nstates = 1 and no exit state
    c = ni(params.cs, .3);
    pientry = ni(params.pi_entries, 3);
    piexit = 1;
%      dl0 = ni(params.dls, 100);
    dl0 = 1;
elseif params.model == "entryexit"
    %good solution for nstates = 5 and entry with exit
     c = ni(params.cs, 10);
    pientry = ni(params.pi_entries, 1);
    piexit = ni(params.pi_exits, .1);
    dl0 = ni(params.dls, 100);
end

%dls, kds, pi_exits, cs, pi_entries
f = factive(dl0:end, kd, piexit, c, pientry);
onset = mfpts(dl0:end, kd, piexit, c, pientry);
dl = params.dls(dl0:end);

f2 = factive(dl0:end, kd2, piexit, c, pientry);
onset2 = mfpts(dl0:end, kd2, piexit, c, pientry); 
f3 = factive(dl0:end, kd3, piexit, c, pientry);
onset3 = mfpts(dl0:end, kd3, piexit, c, pientry); 

figure;
tiledlayout('flow')
nexttile
plot(dl, f)
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
xlim([dl0, 4000])
ylim([0, 1])
nexttile
plot(dl, onset)
xlabel('Dorsal (a.u.)')
ylabel('Mean transcription onset time (min)')
xlim([dl0, 4000])
ylim([0, 8.2])
nexttile
% [X,Y] = meshgrid(dl, f);
% surf(X, Y, onset)
plot3(dl, f, onset, '-o', 'LineWidth', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
grid on
xlim([dl0, 4000])
zlim([0, 8.2])
ylim([0, 1])
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
zlabel('Mean transcription onset time (min)')

% nexttile
figure;
yyaxis left
plot(dl, f, 'LineWidth', 2, 'DisplayName', 'KD=1E3')
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
hold on
plot(dl, f2,  'LineWidth', 2, 'DisplayName', 'KD=1E4')
plot(dl, f3,  'LineWidth', 2,'DisplayName', 'KD=1E5')
xlim([dl0, 4000])
ylim([0, 1])
yyaxis right
plot(dl, onset,  'LineWidth', 2, 'DisplayName', 'KD=1E3')
hold on
plot(dl, onset2,  'LineWidth', 2, 'DisplayName', 'KD=1E4')
plot(dl, onset3,  'LineWidth', 2, 'DisplayName', 'KD=1E5')
xlabel('Dorsal (a.u.)')
ylabel('Mean transcription onset time (min)')
xlim([dl0, 4000])
ylim([0, 8.2])
legend();
title('Delayed kinetic barrier model, one state. c = 2 /min'); 



%%
figure
% for i = 1:size(factive, 2)
for i = kd
    for j = 1:size(factive, 3)
        for k = 1:size(factive, 4)
            for m = 1:size(factive, 5)
                if rand(1) > .5
                    plot3(dl, factive(:, i, j, k, m), mfpts(:, i, j, k, m));
                    hold on
                end
            end
        end
    end
end
grid on
xlim([0, 4000])
zlim([0, 8.2])
ylim([.05, 1])
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
zlabel('Mean transcription onset time (min)')

%%

figure
cmap = viridis(nPlots);
colormap(cmap)
for i = kd
    for j = piexit
        for k = c
            for m = 1:nPlots
                plot3(dl, factive(:, i, j, k, m), mfpts(:, i, j, k, m), 'Color', cmap(m, :), 'LineWidth', 3);
                hold on
            end
        end
    end
end
grid on
xlim([0, 4000])
zlim([0, 8.2])
ylim([0.05, 1])
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
zlabel('Mean transcription onset time (min)')
title('Varying pi_entry')
cb = colorbar;
caxis([min(params.pi_entries), max(params.pi_entries)])
set(gca,'ColorScale','log')
title(cb, 'pi entry')


fmesh = squeeze(factive(:, kd, piexit, c, :));
omesh = squeeze(mfpts(:, kd, piexit, c, :));
dlmesh = repmat(dl', 1, size(fmesh, 2));

hold on

surf(dlmesh, fmesh, omesh, 'FaceAlpha', 0.5, 'FaceColor', [105, 105, 105]/255, 'LineStyle', 'none')



%%
if length(params.pi_exits) > 1
    figure
    cmap = viridis(nPlots);
    colormap(cmap)
    for i = kd
        for j = 1:nPlots
            for k = c
                for m = pientry
                    plot3(dl, factive(:, i, j, k, m), mfpts(:, i, j, k, m), 'Color', cmap(j, :), 'LineWidth', 3 );
                    hold on
                end
            end
        end
    end
    grid on
    xlim([0, 4000])
    zlim([0, 8.2])
    ylim([0.05, 1])
    xlabel('Dorsal (a.u.)')
    ylabel('Fraction of active nuclei')
    zlabel('Mean transcription onset time (min)')
    title('Varying pi_exit')
    cb = colorbar;
    caxis([min(params.pi_exits), max(params.pi_exits)])
    set(gca,'ColorScale','log')
    title(cb, 'pi exit')
    
    
    fmesh = squeeze(factive(:, kd, :, c, pientry));
    omesh = squeeze(mfpts(:, kd, :, c, pientry));
    dlmesh = repmat(dl', 1, size(fmesh, 2));
    
    hold on
    
    surf(dlmesh, fmesh, omesh, 'FaceAlpha', 0.5, 'FaceColor', [105, 105, 105]/255, 'LineStyle', 'none')
end

%%


figure

cmap = viridis(nPlots);
colormap(cmap)

for i = kd
    for j = piexit
        for k = 1:nPlots
            for m = pientry
                plot3(dl, factive(:, i, j, k, m), mfpts(:, i, j, k, m), 'Color', cmap(k, :), 'LineWidth', 3);
                hold on
            end
        end
    end
end
grid on
xlim([0, 4000])
zlim([0, 8.2])
ylim([0.05, 1])
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
zlabel('Mean transcription onset time (min)')
title('Varying c')
cb = colorbar;
caxis([min(params.cs), max(params.cs)])
set(gca,'ColorScale','log')
title(cb, 'c')

fmesh = squeeze(factive(:, kd, piexit, :, pientry));
omesh = squeeze(mfpts(:, kd, piexit, :, pientry));
dlmesh = repmat(dl', 1, size(fmesh, 2));

hold on

surf(dlmesh, fmesh, omesh, 'FaceAlpha', 0.5, 'FaceColor', [105, 105, 105]/255, 'LineStyle', 'none')


%%

figure
% for i = 1:size(factive, 2)
for i = kd
    for j = 1:size(factive, 3)
        for k = 1:size(factive, 4)
            for m = 1:size(factive, 5)
                if rand(1) > .5
                    plot3(dl, factive(:, i, j, k, m), mfpts(:, i, j, k, m));
                    hold on
                end
            end
        end
    end
end
grid on
xlim([0, 4000])
zlim([0, 8.2])
ylim([.05, 1])
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
zlabel('Mean transcription onset time (min)')

%%

% pientry = ni(params.pi_entries, .1);
% c = 10;

figure
cmap = viridis(nPlots);
colormap(cmap)
for i = 1:nPlots
    for j = piexit
        for k = c
            for m = pientry
                plot3(dl, factive(:, i, j, k, m), mfpts(:, i, j, k, m), 'Color', cmap(i, :), 'LineWidth', 3);
                hold on
            end
        end
    end
end
grid on
xlim([0, 4000])
zlim([0, 8.2])
ylim([0.05, 1])
xlabel('Dorsal (a.u.)')
ylabel('Fraction of active nuclei')
zlabel('Mean transcription onset time (min)')
title('Varying KD')
cb = colorbar;
caxis([min(params.kds), max(params.kds)])
set(gca,'ColorScale','log')
title(cb, 'KD')


fmesh = squeeze(factive(:, :, piexit, c, pientry));
omesh = squeeze(mfpts(:, :, piexit, c, pientry));
dlmesh = repmat(dl', 1, size(fmesh, 2));

hold on

surf(dlmesh, fmesh, omesh, 'FaceAlpha', 0.5, 'FaceColor', [105, 105, 105]/255, 'LineStyle', 'none')



%%
% 
% figure;
% t = tiledlayout('flow');
% for k = 1:1:length(params.pi_exits)
%     nexttile;
%     try
%         plotTFDrivenParams2(factive(:, :, :, :, k, :, :),...
%             dt(:, :, :, :, k, :, :), mfpts(:, :, :, :, k, :, :), 'params', params, 'fig', gcf,  'nPoints', 1E3);
%     end
%     %         title(['\pi_{entry} = ', num2str(round2(pi_entries(k))), ' min^{-1}'])
%     title(num2str(round2(params.pi_exits(k))))
% end
% title(t, 'Effect of pi_entry on parameter space');




%% this is for the states heatmap
NumStates = 5;
OnColor = [.6 .96 .64];
Palette = flip((cbrewer('seq', 'Reds',NumStates)));
Palette = ((cbrewer('div', 'Spectral',10)));
Palette = flip((cbrewer('seq', 'Blues',6)));

%Palette = flip(cbrewer('seq', 'YlGnBu',NumStates));
FirstColor = [Palette(1,:) + Palette(2,:)]./2
Palette2 = [Palette([1:6],:);OnColor];
Palette3 = [FirstColor;Palette2([4:6],:);OnColor];
colormap(Palette3)

%% this is for the states heatmap
NumStates = 6;
OnColor = [.6 .96 .64];
Palette = flip(cbrewer('seq', 'OrRd',NumStates));
%Palette = flip(cbrewer('seq', 'YlGnBu',NumStates));
Palette = [Palette;OnColor];
colormap(Palette)

%%
%% this is for the states heatmap
NumStates = 5;
%OnColor = [.6 .96 .64];
Palette = flip((cbrewer('seq', 'Reds',NumStates)));
Palette = ((cbrewer('div', 'RdBu',5)));
Palette = Palette(:,:)
colormap(Palette)




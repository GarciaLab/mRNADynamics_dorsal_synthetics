%%
SilentColor = [0 0 0];
OnColor = [0 1 0];
EntryStatesPalette = flip(cbrewer('seq', 'Greys', 10));
OffStatesPalette = cbrewer('seq', 'OrRd', 5);

FullPalette = [EntryStatesPalette;OffStatesPalette;OnColor;SilentColor];
colormap(FullPalette)

%%
NormalPalette = viridis(12);
NormalPalette = cbrewer('div', 'RdYlBu', 12)
colormap(NormalPalette)
%%
NormalPalette = flip(cbrewer('div', 'RdYlBu', 10))
OnColor = [.6 .9 .6];
SilentColor = NormalPalette(1,:);
FullPalette = [NormalPalette;OnColor;SilentColor];
colormap(FullPalette)

%%
NormalPalette = flip(cbrewer('seq', 'YlGnBu', 12));
SilentColor = [.4 0 .1];
NormalPalette(end,:) = SilentColor;
colormap(NormalPalette)



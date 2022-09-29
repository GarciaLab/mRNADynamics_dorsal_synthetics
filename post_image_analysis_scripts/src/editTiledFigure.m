function editTiledFigure()

fig = gcf
tiles = fig.Children
ax = tiles.Children

for k = 1:length(ax)
    
   set(ax(k), 'YTickLabels', 'hi bud')
   ylabel(ax(k), 'hi bud')
   
end

end
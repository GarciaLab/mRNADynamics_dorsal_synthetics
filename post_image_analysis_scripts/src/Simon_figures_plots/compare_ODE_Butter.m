

%KDs = [100:2500:40000];
KDs= 10:1000:10000;
c=1.5;
dls = linspace(0,4500,20);
Palette = viridis(length(KDs));

figFrac = figure;
axFrac = axes(figFrac);
figOn = figure;
axOn = axes(figOn);

hold(axFrac,'on')
hold(axOn,'on')

for k = 1:length(KDs)
    kd = KDs(k);
    [fraction1,meanOnset1] = fiveSteps(c,dls,kd);
    [fraction2,meanOnset2] = BasicModel_masterEq(c,kd,dls);

    
    plot(axFrac,dls,fraction1,'-o','Color',Palette(k,:))
    plot(axFrac,dls,fraction2,'-*','Color',Palette(k,:))
    ylim([0 1])
    
    plot(axOn,dls,meanOnset1,'-o','Color',Palette(k,:))
    plot(axOn,dls,meanOnset2,'-*','Color',Palette(k,:))
    ylim([0 8])
end


hold(axFrac,'off')
hold(axOn,'off')


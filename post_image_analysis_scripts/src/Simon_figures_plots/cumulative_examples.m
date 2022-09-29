
%% a constant
fig1 = figure;
ax1 = axes(fig1);
fig2 = figure;
ax2 = axes(fig2);


Palette = viridis(7);
Palette = flip(Palette)
Xvals = linspace(0,4000,100);
Constants = logspace(1,1.1,7);

for i = 1:7   
    Yvals = ones(1,length(Xvals)).*Constants(i);
    CumYvals = cumsum(Yvals);
    plot(ax1,Xvals,Yvals,'Color',Palette(i,:),'LineWidth',2)
    hold on
    plot(ax2,Xvals,CumYvals,'Color',Palette(i,:),'LineWidth',2)
    hold on
    
end
hold off

%% hill

figure
Palette = viridis(7);
Palette = flip(Palette)
Xvals = linspace(0,4000,100);
Constants = logspace(3,5,7);
Constants = flip(Constants);
hold on
for i = 1:7
    Kd = Constants(i);
    Yvals = (Xvals/Kd)./(1+(Xvals/Kd));
    Yvals = cumsum(Yvals)
    plot(Xvals,Yvals,'Color',Palette(i,:),'LineWidth',2)
    
end
hold off

%% parabola (repression)

%% hill

figure
Palette = viridis(7);
Palette = flip(Palette)
Xvals = linspace(0,4000,100);
Constants = logspace(3,5,7);
Constants = flip(Constants);
hold on
Constants = linspace(500,1000,7);

for i = 1:7
    Kd = Constants(i);
    Yvals = -0.2.*(Xvals.^2) + Xvals.*Constants(i);
    Yvals = cumsum(Yvals)
    plot(Xvals,Yvals,'Color',Palette(i,:),'LineWidth',2)
    
end
hold off
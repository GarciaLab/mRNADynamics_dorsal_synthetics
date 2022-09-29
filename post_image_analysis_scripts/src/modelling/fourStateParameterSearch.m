R = 200; %aus/min

dl_au = 10:50:3000; %au
co_au = 500; %averageish [Dl]
dl = dl_au/co_au; %dimensionless
pol = 20; %no idea what unit this is

dmrnadts = [];
mrnas = [];

% k01 = 3.77*pol; %pol kon
% k23 = 3.77*pol; %pol on with dl on
k01range = logspace(0, 5);
k23range = logspace(0, 5);
k02range = logspace(0, 5);
k13range = logspace(0, 5);
k10range= logspace(0, 5); %pol off
k32range= logspace(0, 5); %pol off with dl on
k20range= logspace(0, 5); %dl off
k31range= logspace(0, 5); %dl off with pol on

nSamples = 1000;
steepness = nan(1, nSamples);
position = nan(1, nSamples);

for k = 1:nSamples %random samples from parameter space
    
    %         k01= datasample(k01range,1);
    %         k23= datasample(k23range,1);
    %         k02= datasample(k02range,1);
    %         k13= datasample(k13range,1);
    k02coeff = datasample(k02range, 1);
    k13coeff = datasample(k13range, 1);
    k10= datasample(k10range,1);
    k32= datasample(k32range,1);
    k20= datasample(k20range,1);
    k31= datasample(k31range,1);
    
    curve = [];
    temp = [];
    parfor d = 1:length(dl)
        
        k02 = k02coeff*dl(d); %dl on
        k13 = k13coeff*dl(d); %dl on with pol on
        
        
        K = [(-k10-k20) k01 k02 0
            k10 (-k01-k31) 0 k13
            k20 0 (-k02-k32) k23
            0 k31 k32 (-k13-k23)];
        
        
        [t, y, dmrnadt, mrna] = solveFour(K, R, false);
        
        [~, t0] = min(abs(t-6));
        
%         dmrnadts(d, k, :) = dmrnadt;
        temp(d, :) = dmrnadt;
        curve(d) = dmrnadt(t0);
        
    end
    
    dmrnadts(k, :, :) = temp; 
    
    [m, ind] = max(curve);
    curve_norm = curve ./ m;
    [~, dl0] = min(abs(dl - dl(ind)/2));
    dl_norm = dl/dl(dl0);
    slope = diff(curve_norm);
    [steepness(k), s_ind] = max(slope);
    position(k) = dl_norm(s_ind);
%     Ks(:, :, k) = K;
    
end

%
figure;
for k = 1:size(dmrnadts, 2)
    plot(dl_au, dmrnadts(k,:,2));
    waitforbuttonpress
end
xlabel('[Dl] (au)')
ylabel('dmRNA/dt at t = 6 mins (au/min)') 

%get a convex hull of the volume of parameter space, compute its area and
%subdivide it into 9 regions
P = [position',steepness'];
Ps.x = position';
Ps.y = steepness';
[k,av] = convhull(P(:,1),P(:,2));
% [miny, maxyind] = min(steepness)
% [maxy, maxyind] = max(steepness)
% [minx, maxxind] = min(position)
% [maxx, maxxind] = max(position)
% 
% nSlices = 9;
% xinterval = minx:(maxx-minx)/nSlices:maxx;
% yinterval = miny:(maxy-miny)/nSlices:maxy;
% 
% in = [];
% on = [];
% for k = 1:nSlices
%     r{k} = [xinterval(k), yinterval(k); xinterval(k+1), yinterval(k+1)];
%     for j = 1:length(P)
%         [in(j, k),on(j, k)] = inpolygon(P(j, 1),P(j, 2),r{k}(:, 1),r{k}(:,2))
%     end
% end


figure;
tiledlayout('flow')

nexttile;
p = plot(position, steepness, 'o');
hold on
plot(P(k,1),P(k,2), 'r');
fill( position(k),steepness(k), 'g','facealpha', 0.5 );
xlabel('Boundary Position')
ylabel('Steepness')


% 
% nexttile;
% p = plot(position, steepness, 'o');
% hold on
% plot(P(k,1),P(k,2), 'r');
% 
% n = 0;
% for i=1:1:NX
%     for j=1:1:NY
%         n = n + 1;
%         if not(isempty(PXY{i,j}))
%             plot([PXY{i,j}.x PXY{i,j}.x(1)],[PXY{i,j}.y PXY{i,j}.y(1)],'ro-');
%         end
%         hold on
%     end
% end
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% xlabel('Boundary Position')
% ylabel('Steepness')
% % StandardFigurePBoC(p,gca)


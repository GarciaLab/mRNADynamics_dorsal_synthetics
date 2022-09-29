function plotGoodCurves(factive, dt, mfpts, params, varargin)

 %dls, kds, pi_exits, cs, pi_entries

if length(varargin)==1
    goodMatrixIndices = varargin{1};
else
    goodMatrixIndices = plotTFDrivenParams(factive, dt, mfpts, 'dim', 2, 'params', params);
end
 
%sort to get dorsals arrayed contiguously
goodMatrixIndices = sortrows(goodMatrixIndices,[5 4 3 2]);

%get locations of contiguous dorsal arrays
% lenContig = round(length(params.dls)/2); %this needs to be dehardcoded 
% lenContig = length(params.dls);
% lenContig = findContigLength(params, goodMatrixIndices);
lenContig = 8
pos = findArray(goodMatrixIndices(:, 1), lenContig);

figure;
tiledlayout('flow')
 
for j = 1:length(pos)
    
    inds = pos(j) : pos(j)+lenContig-1;
    gmi = goodMatrixIndices( inds, : );
    
    if max(gmi(:, 1)) == length(params.dls) %solns that don't reach max [dl] are weird.

        factive_theory = nan(1, lenContig); 
        dt_theory = nan(1, lenContig); 
        onset_theory = nan(1, lenContig);

        for k = 1:lenContig
            factive_theory(k) = factive(gmi(k, 1), gmi(k, 2), gmi(k, 3), gmi(k, 4), gmi(k, 5));
            dt_theory(k) = dt(gmi(k, 1), gmi(k, 2), gmi(k, 3), gmi(k, 4), gmi(k, 5));
            onset_theory(k) = mfpts(gmi(k, 1), gmi(k, 2), gmi(k, 3), gmi(k, 4), gmi(k, 5));
        end

        if factive_theory(1) < .2 && factive_theory(end) > .6 &&...
               factive_theory(3) < .8
            %only plot curves that span the full factive range
            x = params.dls(gmi(:, 1)); 

            nexttile(1)
            plot(x, factive_theory, 'LineWidth', 2)
            xlim([0, 4000])
            xlabel('[Dorsal] (au)')
            ylabel('fraction of active nuclei')
            ylim([0, 1])
            hold on

            nexttile(2)
            plot(x, dt_theory, 'LineWidth', 2)
            xlabel('[Dorsal] (au)')
            ylabel('Change in mean turn on time across large range of affinities (min)')
            xlim([0, 4000])
            hold on


            nexttile(3)
            plot(x, onset_theory, 'LineWidth', 2)
            hold on
            xlim([0, 4000])
            set(gca, 'XScale', 'log');
            ylim([0, 10])
            ylabel('mean time to turn on (min)')
            xlabel('[Dorsal] (au)')
            
            nexttile(4)
            plot(factive_theory, onset_theory, 'LineWidth', 2)
            hold on
            xlim([0, 1])
            set(gca, 'XScale', 'log');
            ylim([0, 10])
            ylabel('mean transcriptional onset time (min)')
            xlabel('fraction of active nuclei')
            
        end
    end

end


end


function lenContig = findContigLength(params, goodMatrixIndices)

n = zeros(1, length(params.dls));

for k = 1:length(params.dls)
    pos = findArray(goodMatrixIndices(:, 1), k);
    for j = 1:length(pos)
        inds = pos(j) : pos(j)+k-1;
        gmi = goodMatrixIndices( inds, : );
        n(k) = n(k) + (max(gmi(:, 1)) == length(params.dls));
    end
end

[~,lenContig] = min(n(n>length(params.dls)/4));
% lenContig = 22

end
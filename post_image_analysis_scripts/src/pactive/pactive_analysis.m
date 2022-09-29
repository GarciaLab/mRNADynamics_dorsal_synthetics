function pactive_analysis()

load("C:\Users\owner\Dropbox\DorsalSyntheticsDropbox\1Dg11_2xDl\combinedCompiledProjects.mat")

%we're gonna look at what happens to pactive when we remove every other frame from
%traces

prefixes = unique({combinedCompiledProjects.prefix});

p_observeds = [];
p_actives_odd = [];
p_actives_even = [];
Ns = [];

for prefix = 1:length(prefixes)
    
    relation = combinedCompiledProjects( strcmpi({combinedCompiledProjects.prefix}, prefixes{prefix})...
        & [combinedCompiledProjects.cycle] == 12);
    
    N = length(relation);
    Ns(prefix) = N;
    
    oddRel = relation;
    evenRel = relation;
    
    
    
    for k = 1:length(relation)
        evenTemp = [];
        oddTemp = [];
        for j = 1:length(relation(k).particleFrames)
            if mod(relation(k).particleFrames(j), 2)
                oddTemp = [oddTemp, relation(k).particleFrames(j)];
            else
                evenTemp = [evenTemp, relation(k).particleFrames(j)];
            end
        end
        oddRel(k).particleFrames = oddTemp;
        evenRel(k).particleFrames = evenTemp;
    end
    
    
    n_active = 0;
    n_active_odd = 0;
    n_active_even = 0;
    for k = 1:length(relation)
        
        if ~isempty(relation(k).particleFrames)
            n_active = n_active + 1;
        end
        if ~isempty(oddRel(k).particleFrames)
            n_active_odd = n_active_odd + 1;
        end
        if ~isempty(evenRel(k).particleFrames)
            n_active_even = n_active_even + 1;
        end
        
    end
    
    p_active = n_active/N;
    p_active_odd = n_active_odd/N;
    p_active_even = n_active_even/N;
    
    p_observeds(prefix) = p_active;
    p_actives_odd(prefix) = p_active_odd;
    p_actives_even(prefix) = p_active_even;
    
    n_actives(prefix) = n_active;
    n_actives_odd(prefix) = n_active_odd;
    n_actives_even(prefix) = n_active_even;
    
    mean_Dorsals(prefix) = nanmean([relation.dorsalFluoFeature]);
    
end

%I'm only setting these to 0 temporarily because I know for this dataset,
%they are 0. If I use this on other datasets, I need to add a way to keep
%track of these in the loops above and replace these with real values. 
n_actives_evenOnly = 0;
n_actives_oddOnly = 0;


p_detect_even = 1 - (n_actives_evenOnly./n_actives);
p_detect_odd = 1 - ( n_actives_oddOnly./n_actives);

p_detect = (p_detect_even + p_detect_odd) / 2;
p_actives_true = p_observeds ./ p_detect;

c1 = [115,142,193]/255;
c2 = [208,109,171]/255;

plt1 = plot(mean_Dorsals, p_detect, '.');
ax1 = gca;
hold on
plt2 = plot(mean_Dorsals, p_actives_true, '.');
ax2 = gca;
leg = legend('p_detect', 'p_active',  'Interpreter', 'none');
ylim([0, 1.2])
ylabel('fraction')
xlabel('mean Dorsal-Venus fluorescence (au)')
title('Estimated detection probability and estimated fraction competent in 1Dg x 2xDl, nc12')
% StandardFigurePBoC(plt1, ax1)
% StandardFigurePBoC(plt2, ax2)
standardizeFigure(gca, leg);

end
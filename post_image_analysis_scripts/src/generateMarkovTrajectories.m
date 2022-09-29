function generateMarkovTrajectories(whichTransition, tau_entry_off, tau_exit)

%whichTransition is n states x n simulations. a value of 1 corresponds to a
%forward transition from entry -> off -> on. a value of 2 corresponds to a
%transition to the silent state.

nNuclei = min(50, size(whichTransition, 2)); %
nStates = size(whichTransition, 1) + 1;
state_mat = repmat((1:nStates)', [1, nNuclei]); 
time_mat = repmat((zeros(1, nStates))', [1, nNuclei]); 

for k = 1:nNuclei
    
    traj = whichTransition(:, k)';
    transitionToSilent = find(traj==2, 1, 'first');
    
    if ~isempty(transitionToSilent)
        
        state_mat(transitionToSilent+1:end, k) = 12;
        
        time_mat(2:transitionToSilent, k) = cumsum(tau_entry_off(1:transitionToSilent-1, k)); 
        
        time_mat(transitionToSilent + 1:end, k) =...
            ...
            time_mat(transitionToSilent, k) +...
            cumsum(tau_exit(transitionToSilent:end, k));
        
        %let's have it _never_ hop to the off state. once it's on, it stays
        %on
        if state_mat(end-1, k) == nStates-1
            state_mat(end, k) = nStates-1;
        end
        
    else
        state_mat(end, k) = nStates - 1; %let's leave it in the ON state forever (second to last)
        time_mat(2:end, k) = cumsum(tau_entry_off(1:end, k)); 
    end
end


colormap(viridis(size(state_mat', 1)));
figure; imagesc(state_mat')
colorbar;

% T = logspace(-1, log10(8), 100);
T = linspace(0, 8, 10000);

state_matVsTime = ones(nNuclei, length(T));

for n = 1:nNuclei
    for s = 1:nStates
        for t = 1:length(T)-1
            if time_mat(s, n) >= T(t) && time_mat(s, n) < T(t+1)
                state_matVsTime(n, t:end) = state_mat(s, n);
            end
        end
    end
end

figure; 
colormap(viridis(nStates));
imagesc(state_matVsTime)
cb = colorbar;
title(cb, 'state')
ylabel('nucleus')
xlabel('time (min)')
nTicks = 5;
xticks(linspace(1,length(T),nTicks))
xticklabelvec = [T(1:length(T)/(nTicks-1):end-1), max(T)];
xticklabels(string(round(xticklabelvec)))
title('simulation state-space trajectories')

function [fpt_on_observed,factive] = averagePaths_tfdriven(nSims, nSteps, pi0, pi1,pi2, onstate, silentstate, t_cycle, firstoffstate)
%subfunction for tfdrivenentryexit

nSteps = 7;

% fpt_on = [];
fpt_on = nan(1, nSims);
%     fpt_on_observed = [];
%     duration = [];
taus = [exprnd(pi0^-1, [1, nSteps, nSims]) %on
    exprnd(pi1^-1, [1, nSteps, nSims]) %exit
    ];

tau_on = squeeze(taus(1, :, :));
tau_exit = squeeze(taus(2, :, :));

if pi2 ~=0
    tau_entry = exprnd(pi2^-1, [5, nSims]);
end

[~, ind] = min([taus(1,:, :); taus(2,:, :)]);

ind = squeeze(ind);

initStates = 1;
initTimes = 0;


for k = 1:nSims
    
    
        
        [states, times] = makePath(nSteps,onstate, silentstate, firstoffstate, tau_on(:, k), tau_exit(:, k), ind(:, k), initStates, initTimes);
        
        % plot(time, states);
        % xlim([0, 10]);
        %     ton =  times(find(states==onstate, 1 ));
        ton = times(states==onstate);
        %     fpt_on = [fpt_on, ton];
        if ~isempty(ton)
            fpt_on(k) = ton;
        end
        
    
end



fpt_on_observed = fpt_on(fpt_on < t_cycle);

factive = length(fpt_on_observed) / nSims;

% if ~isempty(fpt_on_observed)
%     
%     1
% end

%     duration = t_cycle-fpt_on_observed;

%     histogram(fpt_on_observed, 'Normalization', 'pdf');
%     legend(['<\tau>=', num2str(round2(mean(fpt_on_observed))), ' min']);
%     xlabel('\tau (min)')
%     ylabel('probability density function')


end

function [states, times] = makePath(nSteps, onstate, silentstate,firstoffstate, tau_on, tau_exit, ind, states, times)


n = 6;
for step = 1:nSteps
    n = n + 1;
    if ind(step) == 1 && states(n-1) < onstate
        states(n) = states(n-1) + 1;
        times(n) = times(n-1) + tau_on(step);
    elseif ind(step) == 2 || states(n-1) == onstate
        states(n) = silentstate;
        times(n) = times(n-1) + tau_exit(step);
        return;
    end
end

end


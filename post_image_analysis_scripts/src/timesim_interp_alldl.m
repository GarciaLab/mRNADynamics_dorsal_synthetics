function y = timesim_interp_alldl(x, theta, modelOpts)

%theta- c, kd, n, m, pientry, piexit, t_cycle

if isempty(modelOpts)
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = 1E3;
    modelOpts.modelType = 'entryexit';
end

nSims = modelOpts.nSims;
% exitOnlyDuringOffStates = modelOpts.exitOnlyDuringOffStates;

%mcmcpred feeds in x as a struct. otherwise it's an array
if isstruct(x)
    x = x.ydata(:, 1);
end

dls = x(:); 
n_dls = length(dls);
kd = theta(2);
c = theta(1);
nEntryStates = round(theta(3));
nOffStates = round(theta(4));
pi_entry = theta(5);
pi_exit = theta(6);
t_cycle = theta(7);

onState = nEntryStates + nOffStates + 1;
exitState = onState + 1;

for k = 1:n_dls

    onsets = nan(nSims, 1);
for n = 1:nSims
    
    t = 0;
    state = 1;
    while t < t_cycle
        
%         current_time_index = nearestIndex(tq, t);
        
        if state < nEntryStates + 1
            holding_time = exprnd( pi_entry.^-1);
        else
            
            if pi_exit ~= 0
                
                holding_time = exprnd( (c*(dl./kd) ./ (1 + dl./kd) ).^-1);
                exit_time = exprnd( pi_exit.^-1);
                
                t = t + min(holding_time, exit_time); 
                
                if exit_time < holding_time 
                    state = exitState;
                    break;
                end
                
            else
                holding_time = exprnd( (c*(dl./kd) ./ (1 + dl./kd) ).^-1);
                t = t + holding_time;
                state = state + 1;
            end
            
        end
        
        
        if state == onState
            onsets(n) = t;
            break;
        end
        
    end
    
end

onsets(onsets >= t_cycle) = nan;
fraction_active(k) = numel(onsets(~isnan(onsets)))/nSims;
mean_onset(k) = mean(onsets, 'omitnan');
    

end

y = [fraction_active', mean_onset'];


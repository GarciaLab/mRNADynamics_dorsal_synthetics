function [fraction_active, mean_onset] = timesim_interp(time_vec,dls,varargin)

% numBins is the total number of bins to partition the Dorsal fluorescence;
% bin is which Dorsal fluorescence bin we're simulating

dt = 10/60; %s time resolution, should be the same as in time_vec
c = 1.5;
t_cycle = 8;
%time_vec = linspace(0, t_cycle, ceil(t_cycle/dt))';
kd = 1E3;
%dls =  kd*ones(length(ceil(t_cycle/dt)))';
nSims = 5E3;
nOffStates = 5;
nEntryStates = 0;
pi_entry = 15;
pi_exit = 0;


%%
%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


dls = dls(time_vec < t_cycle);
time_vec = time_vec(time_vec < t_cycle);



q = .1;
tq = min(time_vec):q:max(time_vec);
try
    dq = makima(time_vec, dls, tq)';
catch
    dq = dls(1)*ones(length(tq), 1);
end
% figure; plot(time_vec, dls,'o',tq,dq,'.')

onState = nEntryStates + nOffStates + 1;
exitState = onState + 1;

onsets = nan(nSims, 1);
for n = 1:nSims
    
    t = 0;
    state = 1;
    while t < t_cycle
        
        current_time_index = nearestIndex(tq, t);
        
        if state < nEntryStates + 1
            holding_time = exprnd( pi_entry.^-1);
        else
            
            if pi_exit ~= 0
                
                holding_time = exprnd( (c*(dq(current_time_index)./kd) ./ (1 + dq(current_time_index)./kd) ).^-1);
                exit_time = exprnd( pi_exit.^-1);
                
                t = t + min(holding_time, exit_time); 
                
                if exit_time < holding_time 
                    state = exitState;
                    break;
                end
                
            else
                holding_time = exprnd( (c*(dq(current_time_index)./kd) ./ (1 + dq(current_time_index)./kd) ).^-1);
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
fraction_active = numel(onsets(~isnan(onsets)))/nSims;
mean_onset = mean(onsets, 'omitnan');

%disp('debug stop')
end
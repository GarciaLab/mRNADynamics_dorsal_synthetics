function [fraction_active, mean_onset] = timesim(time_vec,dls,varargin)

% numBins is the total number of bins to partition the Dorsal fluorescence;
% bin is which Dorsal fluorescence bin we're simulating

dt = 10/60; %s time resolution, should be the same as in time_vec
c = 1.5;
t_cycle = 8;
%time_vec = linspace(0, t_cycle, ceil(t_cycle/dt))';
kd = 1E3;
%dls =  kd*ones(length(ceil(t_cycle/dt)))';
nSims = 1E4;
nOffStates = 5;

% grab the Dorsal fluos over time 



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


x = DorsalFluoStruct(1).absoluteTime;
y = DorsalFluoStruct(1).meanDorsalFluo;
tq = min(time_vec):.01:max(time_vec);
dq = makima(x,y,tq);

onState = nOffStates + 1;


onsets = nan(nSims, 1);

for n = 1:nSims
    
    state = ones(length(time_vec), 1);
    
%     tau_offs = exprnd( (c*(dls./kd) ./ (1 + dls./kd) ).^-1, [length(dls), 1]);
       s_offs = poissrnd( (c*(dls./kd) ./ (1 + dls./kd) ), [length(dls), 1]);
   
    %jump_times = [];
    
    for t = 1:length(time_vec)
%         
%         if tau_offs(t) < dt
%             state(t:end) = state(t) + 1;
%         end
        
        state(t+1:end) = state(t) + s_offs(t);

        if state(t) >= onState
           onsets(n) = time_vec(t);
           break;
        end
        
    end
    
    
end


fraction_active = numel(onsets(~isnan(onsets)))/nSims;
mean_onset = mean(onsets, 'omitnan');
%figure

%disp('debug stop')
end
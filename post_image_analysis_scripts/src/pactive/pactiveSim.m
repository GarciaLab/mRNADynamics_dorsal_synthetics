function [fractive, take3PactiveNTrials] = pactiveSim(varargin)

pactive = .6;
pdetect = .5;
N = 20;
nTrials = 5;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

active = zeros(1, N);
trial = zeros(nTrials, N);

rng(1, 'twister')
for k = 1:N
    rActive = rand;
    if rActive < pactive
        active(k) = 1;
    end
    if active(k)        
         trial(:, k) = binornd(1,pdetect, [nTrials, 1]);
    end
    anyTrial(k) = any(trial(:, k));
end 

Nactive = sum(active);
fractive = Nactive./N;

NAnyTrial = sum(anyTrial);
trialSum = sum(trial, 2);
avN_trials = mean(trialSum);
if nTrials == 2
    take3NactiveNTrials = -avN_trials.^2 ./ (NAnyTrial - 2*avN_trials);
else
    syms Nact
    take3NactiveNTrials = vpasolve(1-(NAnyTrial/Nact)==(1-(avN_trials/Nact))^nTrials,Nact, [0, Inf]);
    if length(take3NactiveNTrials) == 2
        take3NactiveNTrials = take3NactiveNTrials(2);
    elseif length(take3NactiveNTrials) > 2
        error('what')
    end
end
% take3NactiveNTrials = -NAnyTrial ./ ((1 - (avN_trials./NAnyTrial)).^nTrials - 1);
take3PactiveNTrials = take3NactiveNTrials ./ N;

if isempty(take3PactiveNTrials)
    take3PactiveNTrials = nan;
end
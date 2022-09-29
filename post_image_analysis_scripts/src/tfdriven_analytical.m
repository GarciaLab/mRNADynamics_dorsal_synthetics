
% model = "basic";
% model = "entry";
% model = "exit";
model = "entryexit";

dmax = 5000;
nPlots = 10;

rng(1)

t_cycle = 10; %min

nSteps = 6;
nSims = 1E3;
nOffStates = 5;
nEntryStates = 5;
firstoffstate = nEntryStates+1;
onstate = nEntryStates + nOffStates+1;
silentstate = onstate+1;
nStates = nEntryStates + nOffStates + 1 + 1;

nEntryTransitions = nEntryStates - 1;
nOffTransitions = nOffStates + 1 - 1; 

occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );

if model == "entry"
    dls = logspace(log10(1), log10(dmax), 100);
    kds = logspace(2, 5, nPlots);
    cs = logspace(1, 3, nPlots);
%     cs = 1;
%     pi1s = logspace(-2, 1, nPlots);
    pi1s = 0;
    pi2s = logspace(-2, 1, nPlots);
%     pi2s = 100;
elseif model == "entryexit"
    nSims = 5E2;
%     dls = linspace(1, dmax, 20);
    dls = logspace(log10(1), log10(dmax), 100);
    kds = logspace(2, 6, nPlots);
    cs = logspace(0, 4, nPlots);
%     cs = 1;
    pi1s = logspace(-2, 1, nPlots);
    pi2s = logspace(-2, 1, nPlots);
%     pi2s = 100; 
elseif model == "basic"
    nSims = 1E3;
%     dls = linspace(1, dmax, 20);
    dls = logspace(log10(1), log10(dmax), 20);
    kds = logspace(2, 6, nPlots*3);
    cs = logspace(0, 4, nPlots*3);
    pi1s = 0;
    pi2s = 1E10;
elseif model == "exit"
    nSims = 1E4;
    dls = logspace(log10(1), log10(dmax), 100);
    kds = logspace(2, 6, nPlots);
    cs = logspace(0, 4, nPlots);
    pi1s = logspace(-2, 1, nPlots);
    pi2s = 1E10;
end

clear params;
params.dls = dls;
params.kds = kds;
params.cs = cs;
params.pi1s = pi1s;
params.pi2s = pi2s;
params.model = model; 

mfpts = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));
factive = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));

%dls, kds, pi1s, cs, pi2s
for n = 1:length(pi2s)
    for m = 1:length(cs)
        for k = 1:length(pi1s)
            for i = 1:length(dls)
                for j = 1:length(kds)
                    pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
                    pi1 = pi1s(k); %min-1
                    pi2 = pi2s(n); %min-1
                    
                    
%                     [fpt_on_observed, factive_temp] = averagePaths_entryexit(nSims, nSteps, pi0, pi1, pi2,onstate, silentstate, t_cycle, firstoffstate);
                    

                    mfpts(i, j, k, m, n) = nEntryTransitions.*(1./(pi2+pi1)) + nOffTransitions.*(1./(pi0 + pi1));

                            
                    factive(i,j,k,m,n) = ((pi2./(pi1+pi2)).^nEntryTransitions ) .*...
                        ((pi0./(pi1+pi0)).^nOffTransitions );
                    
                end
            end
        end
    end
end

plotGoodCurves(factive, dt, mfpts, params)

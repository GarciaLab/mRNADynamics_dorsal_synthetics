function fraction_onset = BasicModel_masterEq(dorsalVals,theta,modelOpts)

%theta = [c, kd, [], NOffStates, [],[], TotalTime]

%Solve the master equation for state of promoter before transcription onset

%% set up problem
% Model parameters:
% c=1;    % transition rate 1/min
% kd = 1E3;
% dls = 1000;

if isempty(modelOpts)
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = 30;%1E3;
    modelOpts.modelType = 'entryexit';
end

%mcmcpred feeds in x as a struct. otherwise it's an array
if isstruct(dorsalVals)
    dorsalVals = dorsalVals.ydata(:, 1);
end

%Simulation paramaters
numCells = modelOpts.nSims;   % the total number of nuclei in the simulation
TotalTime =  theta(7);   % end of the simulation
% dt = 0.1;        %Smaller than all time scales in the system
dt = TotalTime/80;
NOffStates = round(theta(4));   %number of states

c = theta(1);
kd = theta(2);

%Create the matrix to store the results
M(1:TotalTime/dt,1:NOffStates+1) = 0; %initialize to zero everywhenre

%Initial conditions:
M(1,1)=numCells;    %everyone is at state 1 initially

time_vec_2 = 0:dt:TotalTime-dt;

fraction_onset = nan(length(dorsalVals), 2);

forCDF = zeros(size(M,1),length(dorsalVals));
%% Do the calculation
for d = 1:length(dorsalVals)
    dls = dorsalVals(d);
    k = (c*(dls./kd) ./ (1 + dls./kd));
    kdt = k*dt;
    
    for t=2:TotalTime/dt % loop over time steps

        %Calculate the evolution of all boxes minus the ones at the edges
        for s=2:NOffStates % loop over states
            M(t,s) = (1-kdt)*M(t-1,s) + kdt*M(t-1,s-1); %stay + enter - leave
        end

        %Calculate the first box
        M(t,1) = (1-kdt)*M(t-1,1);

        %Calculate the last box
        M(t,NOffStates+1) = M(t-1,NOffStates+1) + kdt*M(t-1,NOffStates);
    end

    fraction_onset(d,1) = M(end,end)/numCells;
    yend = M(:,end); %fraction of nuclei in the last state as a function of time
    forCDF(:,d) = yend;
    fraction_onset(d,2) = sum(diff(yend).*time_vec_2(1:end-1)')/sum(diff(yend)); %expected value
    
    
end


%%
figure(1)
Palette = cbrewer('seq', 'YlGn', length(dorsalVals));
t = time_vec_2;
hold on
for bin = 1:length(dorsalVals)
    binColor = Palette(bin,:);
    CDF = (forCDF(:,bin))./numCells;
    plot(t,CDF,'Color',binColor,'LineWidth',2)
end
hold off
legend(strsplit(num2str(floor(dorsalVals))))
set(gca,'Color',[.9 .9 .87])
ylabel({'cumulative probability'})
%ylabel('cumulative probability')
xlabel('spot turn on time (min)')


figure(2)
Palette = cbrewer('seq', 'YlGn', length(dorsalVals));
t = time_vec_2;
hold on
for bin = 1:length(dorsalVals)
    binColor = Palette(bin,:);
    CDF = (forCDF(:,bin))./numCells;
    diffCDF = diff(CDF);
    plot(t(2:end),diffCDF,'Color',binColor,'LineWidth',2)
end
hold off
legend(strsplit(num2str(floor(dorsalVals))))
set(gca,'Color',[.9 .9 .87])
ylabel({'frequency (across all nuclei)'})
%ylabel('cumulative probability')
xlabel('spot turn on time (min)')


figure(3)
Palette = cbrewer('seq', 'YlGn', length(dorsalVals));
t = time_vec_2;
hold on
for bin = 1:length(dorsalVals)
    binColor = Palette(bin,:);
    CDF = (forCDF(:,bin))./numCells;
    diffCDF = diff(CDF./CDF(end));
    plot(t(2:end),diffCDF,'Color',binColor,'LineWidth',2)
   % plot(,'Color',binColor) %plot mean turn on time.
end
hold off
legend(strsplit(num2str(floor(dorsalVals))))
set(gca,'Color',[.9 .9 .87])
ylabel({'frequency (across active nuclei'})
%ylabel('cumulative probability')
xlabel('spot turn on time (min)')




% %% Make a movie
% MVector=0:NumStates;     %This is the vector of bins for the histogram
% figure(1)
% for t=1:TotalTime/dt
%    bar(MVector,M(t,:))
%    ylim([0,100])
%    drawnow          %Force Matlab to draw the plot 
% end

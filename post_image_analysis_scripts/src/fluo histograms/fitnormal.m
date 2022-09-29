function [mu, sigma] = fitnormal(binedges, vals)

DataX = binedges; %fluorescence bins
DataY = vals; %histogram counts in each bin
% make a gaussian symbolic function to minimize
% this left side of the subtraction is a gaussian evaluated in the domain of the real data
% (DataX), the right hand side is the y values in the data (DataY)
Gaussfun = @(z) 1/(z(2)*sqrt(2*pi)) * exp(-(DataX-z(1)).^2./(2*(z(2)^2)))- DataY; 
% make guesses about gaussian parameters z(1) mu z(2) sigma
z0 = [5.5, .5];
lb = [min(DataX), .01];
ub = [max(DataX), max(DataX)];
Guess = lsqnonlin(Gaussfun,z0,lb,ub); %fit by minizing the subtraction between actual Y data and a gaussian
FitGauss = Gaussfun([Guess(1), Guess(2)]) + DataY;
%function output
FittedMean= Guess(1);
mu = FittedMean;
FittedSD = Guess(2);
sigma = FittedSD;
% plot(DataX,FitGauss,'r','LineWidth',2)
% %hold off
% xlabel('fluorescence (a.u)')
% ylabel('frequency')
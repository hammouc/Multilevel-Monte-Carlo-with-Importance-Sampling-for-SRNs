function [m1,md]= catastrophic_coupling_illustration(h1,exp_number,target,M,tend)

% this code simulates coupled levels of the MLMC explicit TL estimator and plots a
% histogram for the difference data

%% Inputs:
% h1: the time step of the coarest level
% exp_number= the example number (1: decay, 
                                  %2: % (Gene transcription and translation,
                                  %4: the Michaelis-Menten enzyme kinetics )
% M: the number of samples 
% target: the target species
% tend: final time

%% Outputs: 
% m1: mean at the coarse level
% md: mean of the difference

[c,zeta,mu,initial]=example(exp_number);
 S = length(initial);
%actual data for each species at time T.
data = zeros([2*S,M]);


%Simulate the coupled tau-leaping processes.  The 
%output is a vector of length 2*S.  The fist S components are for the 
%tau-leaping process with the final step-size.  The second S components
%are for the tau-leaping process with the cruder step size.

[data(:,1:M)] = coupled_explicit(h1,tend,M,exp_number);

%The difference in the processes.  You could also evaluate f(Z_1) - f(Z_2),
%if you are not simply interested in the abundance of a particular
%species.
md = data(target,1:M)-data(S+target,1:M);

m1 = data(target,1:M);



hist(md)

end %The program
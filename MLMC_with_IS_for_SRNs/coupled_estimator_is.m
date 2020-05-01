function [m1,md,kurt]= coupled_estimator_is(h1,exp_number,target,M,tend,delta)

% this code computes  the mean , the variance of MLMC estimator with importance sampling
% sampling


%% Inputs:
% h1: the time step of the coarest level
% exp_number= the example number (1: decay, 2: Gene transcription and translation, 4: the Michaelis-Menten 
                                                                                     % enzyme kinetics )
% M: the number of samples 
% target: the target species
% tend: final time
% delta:  parameter of the importance sampling algorithm 

%% Outputs: 
% m1: mean at the finest level
% md: mean of the difference
% kurt: kurtosis


tic;
[c,zeta,mu,initial]=example(exp_number);
S = length(initial);
%actual data for each species at time T.
data = zeros([2*S,M]);

%Simulate  paths of the coupled tau-leaping processes.  The 
%output is a vector of length 2*S.  The fist S components are for the 
%tau-leaping process with the final step-size.  The second S components
%are for the tau-leaping process with the cruder step size.

[data(:,1:M),lik,kurt,count] =coupled_explicit_is(h1,tend,M,exp_number,target,delta);

%The difference in the processes.  You could also evaluate f(Z_1) - f(Z_2),
%if you are not simply interested in the abundance of a particular
%species. The difference, md, and the mean, m1, are multiplied by the 
%likelihood factor (lik)
md = (data(target,1:M)-data(S+target,1:M)).*lik;

m1 = data(target,1:M).*lik;





end %The program
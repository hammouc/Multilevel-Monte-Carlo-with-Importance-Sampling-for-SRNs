function [m1,v1,m2,v2,md,vd,CPUTime]= coupled_estimator_is_ML(h1,exp_number,target,M,tend,delta)

% this code computes  the mean , the variance of MLMC estimator,and
% the number of needed samples as well the average cost per level  per sample

%% Inputs:
% h1: the time step of the coarest level
% exp_number= the example number (1: decay, 2: Gene transcription and translation, 4: the Michaelis-Menten 
                            
% M: the number of samples 
% target: the target species
% tend: final time
% delta:  parameter of the importance sampling algorithm 

%% Outputs:
% m1: v1: mean and variance of the fine level
% m2,v2 : mean and variance of the coarsest level
% md,vd : mean and variance of the difference data
% CPUTime: the average cost per level  per sample


display('Coupled Tau-Tau')
h1

[c,zeta,mu,initial]=example(exp_number);
S = length(initial);
%actual data for each species at time T.
data = zeros([2*S,M]);

%For this problem, we only care about the difference of the paths
diff = zeros([1,M]);

bt = clock;
[data(:,1:M),lik,kurt,count] =coupled_explicit_is(h1,tend,M,exp_number,target,delta);
CPUTime = etime(clock,bt);

%The difference in the processes.  You could also evaluate f(Z_1) - f(Z_2),
%if you are not simply interested in the abundance of a particular
%species.
diff(1:M) = (data(target,1:M)-data(S+target,1:M)).*lik;
data(target,1:M) = data(target,1:M).*lik;


m1 = mean(data(target,1:M));
v1 = var(data(target,1:M))/M;
m2 = mean(data(S+target,1:M));
v2 = var(data(S+target,1:M))/M;
md = mean(diff(1:M));
vd = var(diff(1:M));
CPUTime=CPUTime/M;



end %The program
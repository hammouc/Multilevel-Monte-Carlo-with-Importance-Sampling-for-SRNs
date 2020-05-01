function [samples,exp_W,V] = opt(tol,L,L0,M0,exp_number,target,tend,delta)

% this code estimate the variance and average cost per level and then based on these estimates
% computes the number of needed samples  per level for a defined   tolerance (tol)

%% Inputs:
% tol: the defined tolerance
% L: the finest level
% L0: the coarst level
% M0: the initial number of samples 
% exp_number= the example number (1: decay, 2: dimerisation, 3: Lotka)
% target: the target species
% tend: final time
% delta:  parameter of the importance sampling algorithm 

%% Outputs:
% samples: the number of samples needed per level
% exp_W: average cost per level
% V: Variance per level 

Num = L-L0+1;
vals = zeros([Num,2]);

%Generate data for middle steps to use in the optimization problem.
parfor i = 1:(L-L0)
    i
    nt = clock;
    [m1,v1,m2,v2,md,vd,CPUTime]= coupled_estimator_is_ML(tend/2^(L-i),exp_number,target,M0,tend,delta)
    endtime  = etime(clock,nt)
    vals(i,:) = [CPUTime,vd]
end

%Generate data for final tau-leap step to use in the optimization problem.
[m,v,CPUTime] = single_estimator_ML(tend/2^(L0),M0,exp_number,target,tend);
vals(Num,:) = [CPUTime,v]

%%%%%%%%
%%% Optimization problem
%%%%%%%%
exp_W=vals(:,1);
V = vals(:,2);

samples = ceil(2 * sqrt(V./exp_W) * sum(sqrt(V.*exp_W)) / tol^2);
end
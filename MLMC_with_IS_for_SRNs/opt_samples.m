function [ samples] =opt_samples(tol,L,L0,exp_number)
%Uses offline estimates of average cost per level and variances per level for MLMC with IS,
% to compute the optimal number of samples (the offline estimates of costs will
%depend on the machine).  estimates should be computed with opt.m

%-Input: - tol: MLMC tolerance 
%        - L: deepest level of MLMC 
%        - L0: coarsest level of MLMC 
%        - exp_number: example of interest 
%Output:  - samples: optimal number of samples per level to be used for
%MLMC

%% estimated average cost per level  and variances for Example 2 with L0=2 up to L = 10, \detla =0.75
%% computed with M=10^6
if exp_number==2
  
exp_W=[0.053 0.027 0.0138 0.007 0.00376 0.002 0.001 0.00054 0.0003 0.00016 0.000034 ]; %average cost per level
V = [0.000000574 0.0000019 0.000006 0.00002 0.00007 0.00024 0.00078 0.0025 0.008 0.0269 24.07];  %variances per level

% Determing the optimal number of samples per level for MLMC, using equal weights between bias
% and statistical errors
samples = ceil(2 * sqrt(V(11-(L-L0):11)./exp_W(11-(L-L0):11)) * sum(sqrt(V(11-(L-L0):11).*exp_W(11-(L-L0):11))) / tol^2);

%% estimated average cost per level  and variances for Example 4 with L0=0 up to L = 10, \detla =0.75
%% computed with M=10^6
elseif  exp_number==4
exp_W=10^(-3)*[45.6 24 12 6 3 1.53 0.7972 0.373 0.2 0.1  0.0575 0.034 0.006]; %average cost per level
V = [0.00000039 0.00000134 0.0000045 0.000015 0.000051 0.00017 0.00057 0.00193 0.0065 0.022 0.073 0.2363 10];  %variances per level

% Determing the optimal number of samples per level for MLMC, using equal weights between bias
% and statistical errors
samples = ceil(2 * sqrt(V(13-(L-L0):13)./exp_W(13-(L-L0):13)) * sum(sqrt(V(13-(L-L0):13).*exp_W(13-(L-L0):13))) / tol^2);

   
end 



end


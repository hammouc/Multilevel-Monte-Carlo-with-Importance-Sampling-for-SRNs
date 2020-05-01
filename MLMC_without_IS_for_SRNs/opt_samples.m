function [ samples] =opt_samples(tol,L,L0,exp_number)
%Uses offline estimates of average cost per level and variances per level for MLMC without IS,
% to compute the optimal number of samples (the offline estimates of costs will
%depend on the machine). The  estimates should be computed with opt.m
%-Input: - tol: MLMC tolerance 
%        - L: deepest level of MLMC 
%        - L0: coarsest level of MLMC 
%        - exp_number: example of interest 
%Output:  - samples: optimal number of samples per level to be used for MLMC
            
%% estimated average cost per level  and variances for Example 2 with L0=2 up to L = 10,
%% computed with M=10^6
if exp_number==2
    exp_W=[ 0.052 0.027 0.014 0.0073 0.00385 0.00195 0.001 0.00054 0.000295 0.000156 0.00003 ]; %average cost per level
    V = [0.0004 0.0007 0.00118 0.0024 0.004757 0.0096 0.0189 0.036 0.071 0.14 24.067];  %variances per level
    % Determing the optimal number of samples per level for MLMC, using equal weights between bias
% and statistical errors
    samples = ceil(2 * sqrt(V(11-(L-L0):11)./exp_W(11-(L-L0):11)) * sum(sqrt(V(11-(L-L0):11).*exp_W(11-(L-L0):11))) / tol^2);

%% estimated average cost per level  and variances for Example 4 with L0=0 up to L = 10,
%% computed with M=10^6
elseif  exp_number==4
 exp_W=[ 0.045 0.022 0.0117 0.0058 0.003 0.00148 0.00075 0.00037 0.00019 0.0001 0.000055 0.00003 0.0000057];%average cost per level
 V = [ 0.00019 0.00048 0.0008 0.0016 0.0034 0.0066 0.013 0.026 0.053 0.109 0.232 0.575 10 ];  %variances per level
    % Determing the optimal number of samples per level for MLMC, using equal weights between bias
% and statistical errors
 samples = ceil(2 * sqrt(V(13-(L-L0):13)./exp_W(13-(L-L0):13)) * sum(sqrt(V(13-(L-L0):13).*exp_W(13-(L-L0):13))) / tol^2);

 


   
end 

        



end


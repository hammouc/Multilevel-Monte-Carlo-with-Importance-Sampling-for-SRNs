% Example 4 model solved by SSA with M paths
function [average,error,cost] = example4_SSA_MC(M)
params =[0.001 0.005 0.01];
IC=[100 100 0 0];
tspan = [0,1];

val=zeros(1,M);
cost_arr=zeros(1,M);
parfor m=1:M
 tic;
% stochastic models:
[t,y] = SSA(@example3,tspan,IC, params);
cost_arr(m)=toc;
val(m)=y(end,3);
end
average=mean(val);
error=1.96*std(val)/sqrt(M);
cost=mean(cost_arr);


function [S, P, K] = example3(params)
%defines the reactions and the stochastic rate constants for:
S = [1 1 0  0    
     0 0 1  0
     0 0 1  0 ]  ;  % How many of each chemical is used as substractes
P = [0 0 1  0  
     1 1 0  0
     1 0  0 1]  ; % How many does each chemical species is made as products
K  = params';
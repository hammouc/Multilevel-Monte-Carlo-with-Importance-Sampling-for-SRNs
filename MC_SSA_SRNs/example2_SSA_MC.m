% Example 2 model solved by SSA with M paths
function [average,error,cost]=example2_SSA_MC(M)

params =[25 10^3 0.001 0.1 1];
IC=[0 0 0];
tspan = [0,1];

val=zeros(1,M);
cost_arr=zeros(1,M);
parfor m=1:M
 tic;
% stochastic models:
[t,y] = SSA(@example2,tspan,IC, params);
cost_arr(m)=toc;
val(m)=y(end,1);
end
average=mean(val);
error=1.96*std(val)/sqrt(M); % statistical error
cost=mean(cost_arr);

function [S, P, K] = example2(params)
%defines the reactions and the stochastic rate constants for:
S = [0  0  0  
     1  0  0
     0  2  0 
     1  0  0 
     0  1  0 ]  ;  % How many of each chemical is used as substractes
P = [1  0  0  
     1  1  0
     0  0  1
     0  0  0 
     0  0  0 ]  ; % How many does each chemical species is made as products
K  = params';
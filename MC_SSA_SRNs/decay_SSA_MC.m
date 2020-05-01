% decay model solved by SSA with M paths
function [average,error,cost] = decay_SSA_MC(M)
params =[1];
IC=[10];
tspan = [0,1];
val=zeros(1,M);
cost_arr=zeros(1,M);
parfor m=1:M
m
 tic;
% stochastic models:
[t,y] = SSA(@decay,tspan,IC, params);
cost_arr(m)=toc;
val(m)=y(end,1);
end
average=mean(val);
error=1.96*std(val)/sqrt(M); % statistical error
cost=mean(cost_arr);

function [S, P, K] = decay(params)
%defines the reactions and the stochastic rate constants for:
S = [1];  % How many of each chemical is used as substractes
P = [0] ; % How many does each chemical species is made as products
K  = params';
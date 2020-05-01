% Example 2 model solved by SSA, plotting trajectories
function results = example2_SSA_trajectory

% solves dimerisation:
clear all; close all;
params =[25 10^3 0.001 0.1 1];
IC=[0 0 0];
tspan = [0,1];
% figure(1)

% % Create figure
 figure1 = figure;

% % Create axes
 axes1 = axes('Parent',figure1,'FontWeight','bold','FontSize',16);
 box(axes1,'on');
 grid(axes1,'on');
hold(axes1,'all');
for M=1:10
% stochastic models:
[t,y] = SSA(@example2,tspan,IC, params);
y(end,:)
t(end,:)
results.esa = [t,y];
% 
% % % Create multiple lines using matrix input to plot
plot1 = plot(t,y(:,1),'Parent',axes1,'LineWidth',2,'LineStyle','-.');
hold on
% % Create title
%  title({'Trajectories of the spe example by the SSA at final time T=1'},...
%   'FontWeight','bold',...
%   'FontSize',16);
% % 
% % Create legend
 legend(axes1,'show');

end

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
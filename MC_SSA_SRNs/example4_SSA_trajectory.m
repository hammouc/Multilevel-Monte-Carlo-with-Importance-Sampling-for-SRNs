% Example 4 model solved by SSA, plotting trajectories
function results = example4_SSA_trajectory

% solves dimerisation:
clear all; close all;
params =[0.001 0.005 0.01];
IC=[100 100 0 0];
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
[t,y] = SSA(@example3,tspan,IC, params);
y(end,:)
t(end,:)
results.esa = [t,y];
% 
% % % Create multiple lines using matrix input to plot
plot1 = plot(t,y(:,3),'Parent',axes1,'LineWidth',2,'LineStyle','-.');
hold on
% % Create title
%  title({'Trajectories of the spe example by the SSA at final time T=1'},...
%   'FontWeight','bold',...
%   'FontSize',16);
% % 
% % Create legend
 legend(axes1,'show');

end

function [S, P, K] = example3(params)
%defines the reactions and the stochastic rate constants for:
S = [1 1 0  0    
     0 0 1  0
     0 0 1  0 ]  ;  % How many of each chemical is used as substractes
P = [0 0 1  0  
     1 1 0  0
     1 0  0 1]  ; % How many does each chemical species is made as products
K  = params';
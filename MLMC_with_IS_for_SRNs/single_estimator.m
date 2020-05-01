function [m] = single_estimator(h1,M,exp_number,target,tend)
% This code simulate single level MC  with explicit tau leap and computes
% mean, variance, average cost for the estimator

%% inputs
% h1: the time step of the coarest level
% exp_number= the example number (1: decay, 2: Gene transcription and translation,4: the Michaelis-Menten 
                                                                                    %enzyme kinetics )
% M: the number of samples 
% target: the target species
% tend: final time

%% outputs
%m: mean 


%Simulate batch paths of the  tau-leaping processes.  
 
[data(:,1:M)] = single_explicit(h1,tend,M,exp_number);
   
  

m = data(target,1:M);

end %The program
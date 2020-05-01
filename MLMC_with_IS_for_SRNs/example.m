function [c,zeta,mu,initial]=example(exp)
% This function provides the parameters related to  each example 
%% inputs:
 % exp: the index of the exmple
%% outputs
 % c: reactions rates vector
 % mu: the stochiometric matrix 
 % zeta: adjusted stochiometric matrix  to be used in the coupling
 % initial: initial   counting number of each species
 
switch exp
    case 1 % the decay example
   c=1;
   zeta =[[ -1.  -1.  0]
        [ -1.  0.  -1.  ]];
   mu=-1;
   initial=10;
   case 2   %Gene transcription and translation example
  c=[25 10^3 0.001 0.1 1];
 
 zeta =[[ 1.  1.  0. 0. 0.  0.   0.  0.  0. -1.  -1.  0.  0.  0.  0.]
  [ 0.  0.  0.  1.  1.  0.  -2.  -2.  0. 0. 0.  0. -1 -1 0 ]
  [ 0.  0.  0.  0.  0.  0.  1.  1.  0. 0. 0. 0.  0. 0. 0.  ]
  [ 1.  0.  1. 0. 0.  0.   0.  0.  0. -1.  0.  -1.  0.  0.  0.]
  [ 0.  0.  0.  1.  0.  1.  -2.  0.  -2. 0. 0.  0. -1 0 -1 ]
  [ 0.  0.  0.  0.  0.  0.  1.  0.  1. 0. 0. 0.  0. 0. 0.  ]];
   mu=[[1 0 0 -1 0]; [0 1 -2 0 -1]; [0 0 1 0 0]];
  initial=[0 0 0];
   
  case 4 %  the Michaelis-Menten enzyme kinetics  example
  c=[0.001 0.005 0.01];
 
 zeta =[[ -1.  -1.  0. 1. 1.  0.   1.  1.  0.  ]
         [ -1.  -1.  0. 1. 1.  0.   0.  0.  0.]
       [ 1.  1.  0. -1. -1.  0.   -1.  -1.  0.]
         [ 0.  0.  0. 0. 0.  0.   1.  1.  0.]
        [ -1.  0.  -1. 1. 0.  1.   1.  0.  1.  ]
         [ -1.  0.  -1. 1. 0.  1.   0.  0.  0.]
       [ 1.  0.  1. -1. 0.  -1.   -1.  0.  -1.]
         [ 0.  0.  0. 0. 0.  0.   1.  0.  1.]];
  mu=[[-1 1 1]; [-1 1 0];[1  -1 -1];[0  0  1]];
  initial=[ 100 100 0 0];
  
 

end
        
        
end
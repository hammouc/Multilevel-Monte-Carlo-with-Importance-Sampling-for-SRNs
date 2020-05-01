function [r_z]=propensity(z,c,exp)
% this function computes the propensity function (using stochastic mass kinetics principle which corresponds
% to the example given as input
%% inputs:
% z: vector of current state
% c: the reaction rates
% exp: the index of our example
%% Outputs:
% r_z: propensity vector

switch exp
    case 1
    r_z= c*z ;
    
     case 2
    r_z(1)= c(1);
    r_z(2,1)= c(2)*z(1); 
    r_z(3,1)= c(3)*z(2)*(z(2)-1) ; 
    r_z(4,1)= c(4)*z(1) ;
    r_z(5,1)= c(5)*z(2) ;
    case 3
    r_z(1)= c(1)*z(3);
    r_z(2,1)= c(2)*z(1); 
    r_z(3,1)= c(3)*z(3) ; 
    r_z(4,1)= c(4)*z(3) ;
    r_z(5,1)= c(5)*z(2) ;
    r_z(6,1)= c(6)*z(1) *z(2);
     case 4
    r_z(1)= c(1)*z(1) *z(2) ;
    r_z(2,1)= c(2)*z(3) ; 
    r_z(3,1)= c(3)*z(3) ; 
        
end
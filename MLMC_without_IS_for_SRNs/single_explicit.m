function [z] = single_explicit(h1,tend,M,exp_number)

% This code simulates M  paths  of SRNs using single level explicit tau leap  

%% inputs
% h1: the time step of the coarest level
% exp_number= the example number (1: decay,  2: % (Gene transcription and translation, 4: the Michaelis-Menten
                                                                                        %enzyme kinetics )
% M: the number of samples 
% tend: final time

%% Outputs
%z: vector containing the counting number of each species at the final time, tend


[c,zeta,mu,initial]=example(exp_number);
%Number of reaction channels in  network.
R = length(c);
S = length(initial);
%Vector of size (num of species X number of vectors being computed)
z = zeros([S,M]);
 
%A very large maximum number of steps before breaking.
maxi = 10^8;
%initialize.
T = 0;
%Initialize x.  
for j = 1:M
    z(:,j) = initial';
end

% r_z are intensities.  Each column corresponds with an independent vector.  
r_z = zeros([R,M]); 

for i=1:maxi-1,
    
    %If Time is passed, then we should break and take as output the vectors
    %we have.
    if T >= tend
    break
    end
    
    for j = 1:M 
    r_z(:,j)=propensity(z(:,j),c,exp_number);
    end

    %get Poisson random numbers.
    P = poissrnd(r_z*h1);

    %update state vector
    z= z + mu*P;
    z= max( z,0);
    
    %Update time.
    T = T + h1;
   

end % for i=1:maxi-1,



end %The program
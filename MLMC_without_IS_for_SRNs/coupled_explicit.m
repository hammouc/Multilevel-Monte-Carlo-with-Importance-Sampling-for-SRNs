function [x] =coupled_explicit(h1,tend,M,exp_number)

 % This function couple and simulates two consecutive levels of M paths 
 % with coarse time step h1 using explicit tau leap method without IS
 
 %% Inputs:
 % h1: The time step size for the coarse level
 % tend: final time
 % M: The number of needed samples of coupled paths
% exp_number= the example number (1: decay, 2: Gene transcription and translation,4: the Michaelis-Menten 
                                                            %enzyme kinetics )
 
 %% outputs:
 % x: Array(2xM)containing  , x(1:S,:) the sample values for the finer level 
 % at final time and x(1+S:2*S,:) the sample values for the coarse level 
 % at final time
 
[c,zeta,mu,initial]=example(exp_number);

%Number of reaction channels in original network.
J= length(c); 
%Total number of jump directions of the coupled process.
Total = 3*J;

%Number of species in original system.
S = length(initial);

%A very large maximum number of steps before breaking.
maxi = 10^8;

%initialize time.
T = 0;

%Initialize x.  

x=repmat(initial',2,M);


% r_x intensities for the finer level, r_z intensities for the coarse level
r_x = zeros([J,M]);
r_z = zeros([J,M]);

%find intensities for process with coarser step-size.
for j = 1:M
r_z(:,j)=propensity(x(S+1:2*S,j),c,exp_number);
end


lambda = zeros([Total,M]);


for i=1:maxi-1,
    
    %If Time is passed, then we should break and take as output the vectors
    %we have.
    if T >= tend
    break
    end

    %Intensity functions for the different realizations of the tau-leaping 
    %processes with the finer step size.  
    for j = 1:M 
    r_x(:,j)=propensity(x(1:S,j),c,exp_number);
    end


    %Calculuate the intensities for the jump processes of the coupled
    %processes.
    for j = 1:M
    for q = 1:J
    ell = (q-1)*3 + 1;
    lambda(ell,j) = min(r_x(q,j),r_z(q,j));
    lambda(ell + 1,j) = r_x(q,j) - lambda(ell,j);
    lambda(ell + 2,j) = r_z(q,j) - lambda(ell,j);
    end
    end

    %get Poisson random numbers.   
    P = poissrnd(lambda*(h1/2));


    %update state vector
    x(1:S,:)= x(1:S,:) + zeta(1:S,:)*P;
    x(1+S:2*S,:)=  x(1+S:2*S,:) + zeta(1+S:2*S,:)*P;

    x(1:S,:)= max( x(1:S,:),0);
    x(1+S:2*S,:)= max(x(1+S:2*S,:),0);


    %update intensities for process with coarse step-size.   
    if mod(i,2) == 0
    for j = 1:M
    r_z(:,j)=propensity(x(S+1:2*S,j),c,exp_number);
    end
    end

    %Update time.
    T = T + h1/2;
    
    
end % for i=1:maxi-1,




end %The program
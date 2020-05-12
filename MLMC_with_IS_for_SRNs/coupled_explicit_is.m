function [x,lik,kurt,count] =coupled_explicit_is(h1,tend,M,exp_number,target,delta)

 % This function couple and simulates two consecutive levels of M paths 
 % with coarse time step h1 using explicit tau leap method with IS
 % paramterized by the parameter delta
 
 %% Inputs:
 % h1: The time step size for the coarse level
 % tend: final time
 % M: The number of needed samples of coupled paths
% exp_number= the example number (1: decay, 2: Gene transcription and translation,4: the Michaelis-Menten 
                                                            %enzyme kinetics )
 
%delta:  parameter of the importance sampling algoithm 
  
 %% outputs:
 % x: Array(2xM)containing  , x(1:S,:) the sample values for the finer level 
 % at final time and x(1+S:2*S,:) the sample values for the coarse level 
 % at final time
 %lik: likelihood factor accross one path 
 % kurt: kurtosis
 % count: number of times one needs importance sampling



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
count=zeros(1,M);
lik=ones([1,M]);


for i=1:maxi-1
    
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
    % first interval of the coarse level: updating
    if mod(i,2) ~= 0
    b=zeros(J,M);
    %Calculuate the intensities for the jump processes of the coupled
    %processes.
    for j = 1:M
    for q = 1:J
     ell = (q-1)*3 + 1;
     % determining the set J_1 in manuscript by boolean J1_cond
     if (exp_number==1 || exp_number==4)
            J1_cond=0;
     elseif  exp_number==2
            J1_cond=(q~=1) && (q~=4);
      end
    if (r_x(q,j)-r_z(q,j)==0) || ((x(target,j)-x(S+target,j))~=0) || J1_cond 

    if (min(r_x(q,j),r_z(q,j))==r_z(q,j))
        lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
        lambda(ell+1,j) = (r_x(q,j)-r_z(q,j))*(h1/2);
        lambda(ell+2,j) = 0;
    else
        lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
        lambda(ell+1,j) = 0;
        lambda(ell+2,j) =(r_z(q,j)-r_x(q,j))*(h1/2); 
    end
    else
      b(q,j)=1; 
      count(j)=count(j)+1;
      if (min(r_x(q,j),r_z(q,j))==r_z(q,j))
          lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
          lambda(ell+1,j) = (r_x(q,j)-r_z(q,j))*(h1/2)^(1-delta);
          lambda(ell+2,j) = 0;
      else
          lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
          lambda(ell+1,j) = 0;
          lambda(ell+2,j) =(r_z(q,j)-r_x(q,j))*(h1/2)^(1-delta); 
      end 

    end
    end
    end

    %get Poisson random numbers.   
    P = poissrnd(lambda);


    %update state vector
    x(1:S,:)= x(1:S,:) + zeta(1:S,:)*P;
    x(1+S:2*S,:)=  x(1+S:2*S,:) + zeta(1+S:2*S,:)*P;

    x(1:S,:)= max( x(1:S,:),0);
    x(1+S:2*S,:)= max(x(1+S:2*S,:),0);

    %Computing the likelihood
    for j = 1:M
    for q=1:J
    ell = (q-1)*3 + 1;
    if b(q,j)==1

        lik(j)= lik(j)*exp((h1/2)*((h1/2)^(-delta)*abs(r_x(q,j)-r_z(q,j))-abs(r_x(q,j)-r_z(q,j))))*((h1/2)^(delta))^(abs(P(ell+1,j)-P(ell+2,j)));
    end
    end
    end 


    % Second interval of the coarse level: updating
    elseif mod(i,2) == 0
    b=zeros(J,M);
    %Calculuate the intensities for the jump processes of the coupled
    %processes.
    for j = 1:M
    for q = 1:J
     ell = (q-1)*3 + 1;
      % determining the set J_1 in manuscript by boolean J1_cond
     if (exp_number==1 || exp_number==4)
            J1_cond=0;
     elseif  exp_number==2
            J1_cond=(q~=1) && (q~=4);
      end
    if (r_x(q,j)-r_z(q,j)==0) || ((x(target,j)-x(S+target,j))~=0) || J1_cond 


    if (min(r_x(q,j),r_z(q,j))==r_z(q,j))
        lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
        lambda(ell+1,j) = (r_x(q,j)-r_z(q,j))*(h1/2);
        lambda(ell+2,j) = 0;
    else
        lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
        lambda(ell+1,j) = 0;
        lambda(ell+2,j) =(r_z(q,j)-r_x(q,j))*(h1/2); 
    end
    else
      b(q,j)=1; 
      count(j)=count(j)+1;
      if (min(r_x(q,j),r_z(q,j))==r_z(q,j))
          lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
          lambda(ell+1,j) = (r_x(q,j)-r_z(q,j))*(h1/2)^(1-delta);
          lambda(ell+2,j) = 0;
      else
          lambda(ell,j) = min(r_x(q,j),r_z(q,j))*(h1/2);
          lambda(ell+1,j) = 0;
          lambda(ell+2,j) =(r_z(q,j)-r_x(q,j))*(h1/2)^(1-delta); 
      end 

    end
    end
    end

    %get Poisson random numbers.   
    P = poissrnd(lambda);


    %update state vector
    x(1:S,:)= x(1:S,:) + zeta(1:S,:)*P;
    x(1+S:2*S,:)=  x(1+S:2*S,:) + zeta(1+S:2*S,:)*P;

    x(1:S,:)= max( x(1:S,:),0);
    x(1+S:2*S,:)= max(x(1+S:2*S,:),0);

    %Computing the likelihood
    for j = 1:M
    for q=1:J
    ell = (q-1)*3 + 1;
    if b(q,j)==1

        lik(j)= lik(j)*exp((h1/2)*((h1/2)^(-delta)*abs(r_x(q,j)-r_z(q,j))-abs(r_x(q,j)-r_z(q,j))))*((h1/2)^(delta))^(abs(P(ell+1,j)-P(ell+2,j)));
    end
    end
    end 

    %update intensities for process with coarse step-size.
    for j = 1:M
    r_z(:,j)=propensity(x(S+1:2*S,j),c,exp_number);
    end

    end

    %Update time.
    T = T + h1/2;
    
    
end % for i=1:maxi-1,



clear lambda;
clear P;
clear r_z
%computing the kurtosis of the difference
kurt=kurtosis((x(target,:)-x(target+S,:)).*lik);


end %The program
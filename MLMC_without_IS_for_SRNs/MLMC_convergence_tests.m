%% MLMC convergence tests 
% This code performs MLMC convergence tests (without importance sampling): saves in .txt file and plots variances, expectations, kurtosis per
% level for MLMC

function MLMC_convergence_tests

close all; clear all;

addpath('..');

N0   = 10000;    % initial samples on coarse levels
Lmin = 2;      % minimum refinement level
Lmax = 10;     % maximum refinement level

% Request the user to enter the example of interest and return an error if
% the example in not in example.m

example = 'What is the example of interest? ';
exp=input(example);
if ~ismember(exp, [1 2 4]) 
 error('Error: Example must be defined in example.m')
end 


 
%% defining the total number of samples as well the deepest level for MLMC, for each example in example.m  
  if (exp==1) 
    fprintf(1,'\n ---- decay_example_explicit_tau_leap---- \n');
    N      = 1000000;    % samples used for convergence tests
    L      = 12;        % total number of levels for convergence tests 
    Eps    = [  0.01 0.02 0.05 0.1];
  elseif (exp==2) 
    fprintf(1,'\n ---- example2_explicit_tau_leap---- \n');
    N      = 1000000;    % samples used for convergence tests
    L      = 10;       % total number of levels for convergence tests 
    Eps    = [ 0.01 0.02 0.05 0.1 ];
  elseif (exp==4) 
    fprintf(1,'\n ---- example4_explicit_tau_leap---- \n');
    N      = 1000000;    % samples used for convergence tests
    L      = 10;         % total number of levels for convergence tests
    Eps    = [ 0.01 0.02 0.05 0.1 ];
  end
 
 %% writing the outputs in .txt file
 
  filename = ['Example_' num2str(exp)];
  fp = fopen([filename '.txt'],'w');
  mlmc_test(@mlmc_tau_leap_l, N,L, N0,Eps,Lmin,Lmax, fp, exp);
  fclose(fp);


%% print  reference value

  if (exp==1)
   val=3.6787944117144235;   % computed analytically
  elseif (exp==2)
   val=23.79;                %% computed with an exact scheme SSA with statistical error ± 0.004
  elseif (exp==4)
   val=9.04;                 %% computed with an exact scheme SSA with statistical error ± 0.004
  end

  if isnan(val)
    fprintf(1,'\n Exact value unknown \n\n');
  else
    fprintf(1,'\n Exact value: %f \n\n',val);
  end

%
%% plot results
%
  nvert = 2;
  mlmc_plot(filename, nvert);

  if(nvert==1)
    figure(1)
    print('-deps2',[filename 'a.eps'])
    figure(2)
    print('-deps2',[filename 'b.eps'])
  else
    print('-deps2',[filename '.eps'])
  end

%
%% now do 100 MLMC calcs in parallel
%
  filename = ['mlmc_tau_leap_l_' num2str(exp) '_100'];
  fp = fopen([filename '.txt'],'w');
  mlmc_test_100(@mlmc_tau_leap_l, val, N0,Eps,Lmin,Lmax, fp, exp);
  fclose(fp);

%
%% plot results
%
  mlmc_plot_100(filename);
  print('-deps2',[filename '.eps'])


%-------------------------------------------------------
%
%% defing the function mlmc_tau_leap_l for the single and coupled levels for MLMC estimator
%  -Inputs: - l: level
%           - N: number of samples used to estimate the moments 
%           - exp: example of interest
%  -Outputs: - sums: related to different moments used for the estimation of
%            convergence rates for weak and strong error and computing the
%             kurtosis
   
function [sums, cost] = mlmc_tau_leap_l(l,N, exp)

% defining the final time (tend), Coarsest level, hierarchy of number of time steps and
% step sizes for coarse and fine level, for each example
if exp==1
    tend=1.0; % defining final time
    nf = 2^(l); % hierarchy of number of time steps for fine levels
    nc = nf/2;   % hierarchy of number of time steps for coarse levels

    hf = tend/nf; % hierarchy of step sizes for fine levels
    hc = tend/nc; % hierarchy of step sizes for coarse levels
elseif  exp==2
    tend=1.0;
    nf = 2^(l+2);
    nc = nf/2;

    hf = tend/nf;
    hc = tend/nc;    
elseif  exp==4
    tend=1.0;
    nf = 2^(l);
    nc = nf/2;

    hf = tend/nf;
    hc = tend/nc;
end
sums(1:6) = 0;

for N1 = 1:10000:N
  N2 = min(10000,N-N1+1);
 
  if exp==1
    if l==0
    Pf = single_estimator(hf,N2,1,1,tend); % single level estimator
    elseif l>0
    [Pf,dP]= coupled_estimator(hc,1,1,N2,tend); % coupled level estimator
    end
  
  elseif  exp==2
    if l==0
    Pf = single_estimator(hf,N2,exp,1,tend);
    elseif l>0
    [Pf,dP]= coupled_estimator(hc,exp,1,N2,tend);
   end  
       
  elseif  exp==4
    if l==0
    Pf = single_estimator(hf,N2,exp,3,tend);
    elseif l>0
    [Pf,dP]= coupled_estimator(hc,exp,3,N2,tend);
    end
  end
        

  if l==0
    dP = Pf;
  end

  sums(1) = sums(1) + sum(dP);
  sums(2) = sums(2) + sum(dP.^2);
  sums(3) = sums(3) + sum(dP.^3);
  sums(4) = sums(4) + sum(dP.^4);
  sums(5) = sums(5) + sum(Pf);
  sums(6) = sums(6) + sum(Pf.^2);
end

cost = N*nf;   % cost defined as number of fine timesteps



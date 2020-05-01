function [cputime,Mean,Variance,samples,cost] = runMLMC_SRN(tol,L,L0,exp_number,target,tend,delta)

% this code computes  the mean , the variance of MLMC estimator with importance sampling, 
%the number of needed samples as well the average cost per level for a defined tolerance (tol)

%% Inputs:
% tol: the defined tolerance
% L: the finest level (you start with a small guess and the algorithm
% increases L until reaching the bias (weak error) constraint
% L0: the coarst level
% M0: the initial number of samples 
% exp_number= the example number (1: decay, 2: dimerisation, 3: Lotka)
% target: the target species
% tend: final time
% delta:  parameter of the importance sampling algorithm 

%% Outputs:
% cputime: the total time of the MLMC estimator code
% Mean: the mean of MLMC estimator for the defined tolerance Tol
% Variance: the variance of MLMC estimator for the defined Tol
% samples: the number of samples needed per level
% cost: the average cost per level

bt_total = clock; %current date vector
Num = L-L0+1;


samples =opt_samples(tol,L,L0,exp_number);  % use saved estimated of variance and cost 
% per level to determine the optimal number of samples, or use opt below to have a new estimates 
%alternatively
%M0=10^6;
%samples= opt(tol,L,L0,M0,exp_number,target,tend,delta);


%%%%%%%%
%Now to run the actual MLMC.
%%%%%%%%

%mean_diff will track the differences in the mean (for the coupled processes), 
%and the mean of the tau-leap process (coursest level).  variances will track
%the variances at each level. 
mean_diff = zeros([Num,1]);
variances = zeros([Num,1]);
cost= zeros([Num,1]);
%clock to see how long this inner MLMC takes.  This is the bulk of the
%code and where most work will take place.
bt_inner = clock;
%Generate data for middle steps until reaching bias (weak error) constraint
convergence=0;
while (convergence==0)
  
     sample_new = samples(1);
   [m1,v1,m2,v2,md,vd,CPUTime]= coupled_estimator_is_ML(tend/2^(L-1),exp_number,target,sample_new,tend,delta);
     mean_diff(1) = md;
     variances(1) = vd/sample_new;
     md1=md;
     cost(1)=CPUTime;
     
      sample_new = samples(2);
    [m1,v1,m2,v2,md,vd,CPUTime]= coupled_estimator_is_ML(tend/2^(L-2),exp_number,target,sample_new,tend,delta);
     mean_diff(2) = md;
     variances(2) = vd/sample_new;
     cost(2)=CPUTime;
     md2=md;
     if     max(abs(md2)/2,abs(md1))<tol/sqrt(2)
         convergence=1
     else
         L=L+1
        
         samples =opt_samples(tol,L,L0,exp_number); % use saved estimated of variance and cost per 
                        %level to determine the optimal number of samples, or use opt below to have
                        %a new estimates 
        %samples= opt(tol,L,L0,M0,exp_number,target,tend,delta);
         Num = L-L0+1;
     
     end
     
end
for i = 3:(L-L0)

      sample_new = samples(i);

    [m1,v1,m2,v2,md,vd,CPUTime]= coupled_estimator_is_ML(tend/2^(L-i),exp_number,target,sample_new,tend,delta);
    
     mean_diff(i) = md;
     variances(i) = vd/sample_new;
     cost(i)=CPUTime;
    
end

%Generate data for final step.
sample_new = samples(Num);
[m,v,CPUTime] = single_estimator_ML(tend/2^(L0),sample_new,exp_number,target,tend);
mean_diff(Num) = m;
variances(Num) = v/sample_new;
cost(Num)=CPUTime;

%output differences in mean and the variances.
mean_diff;
variances;

Mean = sum(mean_diff);
Variance = sum(variances);

cputime_inner = etime(clock,bt_inner)

cputime = etime(clock,bt_total);

   
end


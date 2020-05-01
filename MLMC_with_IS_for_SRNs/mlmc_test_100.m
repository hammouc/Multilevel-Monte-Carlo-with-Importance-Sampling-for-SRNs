%
% function mlmc_test_100(mlmc_l, val, N0,Eps,Lmin,Lmax, fp, varargin)
%  test routine to perform 100 independent MLMC calculations in parallel
%
% sums = mlmc_l(l,N, varargin)     low-level routine
%
% inputs:  l = level
%          N = number of paths
%
% output: sums(1) = sum(Pf-Pc)
%         sums(2) = sum((Pf-Pc).^2)
%
% val      = exact value (NaN if not known)
% N0       = initial number of samples for MLMC calcs
% Eps      = desired accuracy array for MLMC calcs
% Lmin     = minimum number of levels for MLMC calcs
% Lmax     = maximum number of levels for MLMC calcs
% fp       = file handle for printing to file
% varargin = optional additional user variables to be passed to mlmc_l
%

function mlmc_test_100(mlmc_l, val, N0,Eps,Lmin,Lmax, fp, varargin)

PRINTF2(fp,'\n');
PRINTF2(fp,'**********************************************************\n');
PRINTF2(fp,'*** MLMC_100 file version 0.9     produced by          ***\n');
PRINTF2(fp,'*** MATLAB mlmc_test_100 on %s       ***\n',datestr(now));
PRINTF2(fp,'**********************************************************\n');
PRINTF2(fp,'\n');
PRINTF2(fp,'***************************************** \n');
PRINTF2(fp,'*** MLMC errors from 100 calculations *** \n');
PRINTF2(fp,'***************************************** \n');

if isnan(val)
  PRINTF2(fp,'\n Exact value unknown \n');
else
  PRINTF2(fp,'\n Exact value: %f \n',val);
end

for i = 1:length(Eps)
  PRINTF2(fp,'\n eps = %.3e \n-----------------\n',Eps(i)); 

  P100 = zeros(1,100);
  parfor j=1:100
    RandStream.setGlobalStream( ...
    RandStream.create('mrg32k3a','NumStreams',100,'StreamIndices',j));

    P100(j) = mlmc(mlmc_l,N0,Eps(i),Lmin,Lmax, 0,0,0, varargin{:});
  end
  for j=1:5:100
    PRINTF2(fp,' %.5e ',P100(j:j+4));
    PRINTF2(fp,'\n');
  end
end

end

%
% function to print to both a file and stdout 
%

function PRINTF2(fp,varargin)
fprintf(fp,varargin{:});
fprintf( 1,varargin{:});
end
%
% utility to generate MLMC plots based on input text file  mlmc_plot(filename,nvert)
%

function mlmc_plot(filename,nvert,varargin)

close all;

if nargin==2
  error_bars = 0;
elseif nargin==3
  if varargin{1}=='error_bars'
    error_bars = 1;
  else
    error('invalid mlmc_plot option') 
  end
else
  nargin
  error('invalid number of mlmc_plot arguments') 
end

%
% read in data
%

fid = fopen([filename '.txt'],'r');

line = '    ';
while (length(line)<20) | (strcmp(line(1:4),'*** ')==0)
  line = [ fgetl(fid) '    ' ];
end
file_version = sscanf(line(23:30),'%f');
if isempty(file_version)
  file_version = 0.8;
end

if (file_version<0.9)
  N = 0;
  if (error_bars)
    error('cannot plot error bars -- no value of N in file');
  end
else
  while (length(line)<20) | (strcmp(line(1:9),'*** using')==0)
    line = [ fgetl(fid) '    ' ];
  end
  N = sscanf(line(14:20),'%d');
end

line = '    ';
while (length(line)<20) | (strcmp(line(1),'-')==0)
  line = [ fgetl(fid) '    ' ];
end

line = fgetl(fid);
l    = 1;
while (length(line)>10)
  data = sscanf(line,'%f');
  del1(l) = data(2);
  del2(l) = data(3);
  var1(l) = data(4);
  var2(l) = data(5);
  kur1(l) = data(6);
  chk1(l) = data(7);
  cost(l) = data(8);

  line = fgetl(fid);
  l    = l+1;
end

vvr1 = var1.^2 .* (kur1-1);

L = l-2;

line = '    ';
while (length(line)<20) | (strcmp(line(1),'-')==0)
  line = [ fgetl(fid) '    ' ];
end

line = fgetl(fid);
l    = 1;

while (length(line)>10)
  data = sscanf(line,'%f');
  Eps(l)       = data(1);
  mlmc_cost(l) = data(3);
  std_cost(l)  = data(4);
  len          = length(data)-5;
  ls(1:len,l)  = 0:len-1;
  Nls(1:len,l) = data(6:end);

  line = fgetl(fid);
  l    = l+1;
end

%
% plot figures
%
L1 = 1:L;

figs(1) = figure; 
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75*nvert]; set(gcf,'pos',pos);

set(0,'DefaultAxesColorOrder',[0 0 0]);
set(0,'DefaultAxesLineStyleOrder','-*|:*|--o|--x|--d|--s')
%set(0,'DefaultAxesLineStyleOrder','-*|--*')

subplot(nvert,2,1)
plot(0:L,log2(var2),1:L,log2(var1(2:end)),L1,log2(2.^(-5*L1/4)))
xlabel('level $\ell$','Interpreter','latex'); 
ylabel('$\log_2$ variance','Interpreter','latex');
current = axis; axis([ 0 L current(3:4) ]);
legend({'$g_\ell$','$g_\ell\!-\! g_{\ell-1}$','$2^{-5\ell/4}$'}, ...
        'Interpreter','latex','Location','SouthWest')
hold on;

if error_bars
plot([1:L; 1:L],[log2(max(abs(var1(2:end))-3*sqrt(vvr1(2:end)/N),1e-10)); ...
                 log2(    abs(var1(2:end))+3*sqrt(vvr1(2:end)/N))],'-r.') 
end

subplot(nvert,2,2)
plot(0:L,log2(abs(del2)),1:L,log2(abs(del1(2:end))),L1,log2(2.^(-L1)))
xlabel('level $\ell$','Interpreter','latex'); 
ylabel('$\log_2 |\mbox{mean}|$','Interpreter','latex');
current = axis; axis([ 0 L current(3:4) ]);
legend({'$g_\ell$','$g_\ell\!-\! g_{\ell-1}$','$2^{-\ell}$'}, ...
        'Interpreter','latex','Location','SouthWest')
hold on;

%plot([0:L; 0:L],[log2(max(abs(del2)-3*sqrt(var2/N),1e-10)); ...
%                 log2(    abs(del2)+3*sqrt(var2/N))],'-r.') 
if error_bars
plot([1:L; 1:L],[log2(max(abs(del1(2:end))-3*sqrt(var1(2:end)/N),1e-10)); ...
                 log2(    abs(del1(2:end))+3*sqrt(var1(2:end)/N))],'-r.') 
end

if nvert==2
  subplot(2,2,3)
  %  plot(1:L,chk1(2:end),'--*')
  %  xlabel('level l'); ylabel('consistency check');
  plot(0:L,log2(cost(1:end)),'--*')
  hold on
  plot(L1,log2(2.^(L1)),'--o')
  
  xlabel('level $\ell$','Interpreter','latex'); 
  ylabel('$\log_2$ cost per level','Interpreter','latex');
  legend({'$\log_2$ cost per level','$2^{\ell}$'}, ...
        'Interpreter','latex','Location','NorthWest')
  current = axis; axis([ 0 L current(3:4) ]);

  subplot(2,2,4)
  plot(1:L,kur1(2:end),'--*')
  xlabel('level $\ell$','Interpreter','latex'); 
  ylabel('kurtosis');
  current = axis; axis([ 0 L current(3:4) ]);
end

if nvert==1
  figs(2) = figure;
  pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75]; set(gcf,'pos',pos);
end


set(0,'DefaultAxesLineStyleOrder',':o|:x|:d|:*|:s');
%set(0,'DefaultAxesLineStyleOrder','--o|--x|--d|--*|--s');

if nvert==3
subplot(nvert,2,2*nvert-1)
semilogy(ls, Nls)
xlabel('level $\ell$','Interpreter','latex'); 
ylabel('$N_\ell$','Interpreter','latex'); 
current = axis; axis([ 0 size(Nls,1)-1 current(3:4) ]);
for i=1:length(Eps)
  labels{i} = num2str(Eps(i));
end
legend(labels,'Location','NorthEast')

set(0,'DefaultAxesLineStyleOrder','-*|:*')
%set(0,'DefaultAxesLineStyleOrder','-*|--*')

subplot(nvert,2,2*nvert)
loglog(Eps,Eps.^2.*std_cost(:)', Eps,Eps.^2.*mlmc_cost(:)')
xlabel('accuracy $\varepsilon$','Interpreter','latex'); 
ylabel('$\varepsilon^2$ Cost','Interpreter','latex');
current = axis; axis([ Eps(1) Eps(end) current(3:4) ]);
legend('Std MC','MLMC')
end
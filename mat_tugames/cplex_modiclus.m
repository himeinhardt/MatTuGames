function [x1, fmin]=cplex_modiclus(v,tol)
% CPLEX_MODICLUS computes the modiclus of game v using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
% 
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%
% Usage: [x, alp]=cplex_modiclus(v,tol)
% Define variables:
%  output:
%  x1        -- The modiclus of game v.
%  fmin      -- The minmax bi-excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/19/2017        0.9             hme
%   04/04/2020        1.9             hme
%                


warning('off','all');
if nargin<2
 tol=10^6*eps;
end
%tol=-tol;

N=length(v);
[~, n]=log2(N);

% solver parameter
%ub=[];
%lb=[];
x0=[];
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else
%  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
  options.largescale='on';
  options.algorithm='dual-simplex';
  options.tolfun=1e-10;
  options.tolx=1e-10;
  options.tolrlpfun=1e-10;
  %%%% for dual-simplex
  % opts.MaxTime=9000;
  options.preprocess='none';
  options.tolcon=1e-6;
  options.maxiter=10*(N+n);
  options.display='off';
  options.threads=3;
end
ra = reasonable_outcome(v);
ub=[ra,inf];
lb=-inf(1,n+1);

S=0:N;
N1=N+1;
n1=2*n;
N2=(N1)^2-N1;
for k=1:n, A1(:,k) = -bitget(S,k);end
ve=[0,v];
B1=zeros(N2,1);
A2=zeros(N2,n);  %% Could produce an out of memory!!!
cl=zeros(1,N2);
ii=1;
for k=1:N1
    for jj =1:N1
        if k ~= jj
           if k>1 && jj >1
              cl(ii)=(k-1)+(jj-1)*N1;
           elseif k==1 && jj >1
              cl(ii)=N1*(jj-1);
           elseif k>1 && jj==1
              cl(ii)=k-1;
           end
           A2(ii,:) = A1(k,:)-A1(jj,:);
           B1(ii) = ve(jj)-ve(k);
           ii = ii+1;
        end
    end
end
Aeq=[ones(1,n),0];
A2(:,end+1)=-1;
Beq=v(N);
C=[zeros(n,1);1];
it=0:-1:1-n1;
bA=find(A2(:,end)==0)';
while 1
  A2=sparse(A2);
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,Aeq,Beq,lb,ub,x0,options);
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(lambda.ineqlin>tol))';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     warning('on','all');
     break;
  end
  bA=[bA,bS2];
  S2=cl(bA);
  mS2=rem(floor(S2(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+fmin;
  if rk==n1 
     x=(-mS2\B1(bA))';
     break;
  end
  A2(bS2,end)=0;
end

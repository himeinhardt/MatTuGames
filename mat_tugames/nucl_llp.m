function [x1, fmin]=nucl_llp(v,tol)
% NUCL_LLP computes the nucleolus of game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=nucl_lpp(v,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game v.
%  fmin      -- The minmax excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/23/2014        0.6             hme
%   03/28/2015        0.7             hme
%   06/22/2023        1.9.1           hme
%   04/22/2024        1.9.2           hme
%                


if nargin<2
 tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
if N==3
  x1=StandardSolution(v);
  fmin=-inf;
  return
end

S=1:N;
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[-v';v(N)];
C=[zeros(1,n),1];

k=1:n;
ra = reasonable_outcome(v);
vi=v(bitset(0,k));
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if sum(vi)>v(N)
   error('sum of lower bound exceeds value of grand coalition! No solution can be found that satisfies the constraints.')
end
lb=[vi,-Inf];
ub=[ra,Inf];

options.Display='off';
options.Simplex='on';
options.LargeScale='on';
mth1=verLessThan('matlab','24.1.0');
if mth1==0,
    opts.Algorithm='dual-simplex-highs';
else
    opts.Algorithm='dual-simplex';
end
options.TolFun=1e-10;
options.TolX=1e-10;
options.TolRLPFun=1e-10;
options.MaxIter=256;
%opts
%% for dual-simplex
options.MaxTime=9000;
options.Preprocess='none';
options.TolCon=1e-6;
options.MaxIter=10*(N+n);
it=0:-1:1-n;
bA=find(A1(:,end)==0)';
while 1
  try
    [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,options);
  catch %% old api (before R2022a) with initial value.
   opts.Algorithm='dual-simplex'; 	  
   [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],options);
  end
  x=xmin;
  x1=x';
  x1(end)=[];
  lambda.ineqlin;
  bS1=find(lambda.ineqlin'>tol);
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  bA=[bA,bS2];
  mS2=rem(floor(bA(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+fmin;
  if rk==n 
     x=(-mS2\B1(bA))';
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
end

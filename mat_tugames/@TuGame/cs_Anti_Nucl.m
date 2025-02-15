function [x1, fmin]=cs_Anti_Nucl(clv,cs,tol)
% CS_ANTI_NUCL computes the anti nucleolus of game v w.r.t. coalition structure cs 
% using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=clv.cs_Anti_Nucl(cs,tol)
% Define variables:
%  output:
%  x1        -- The anti nucleolus of game v w.r.t. coalition structure cs. 
%  fmin      -- The maxmin excess value.
%
%  input:
%  clv      -- TuGame class object.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/22/2017        0.9             hme
%   05/26/2024        1.9.2           hme
%                


if nargin<3
 tol=10^6*eps;
end


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

k=1:n;
ra = clv.smallest_amount(); %% Cannot be used as an lower bound. Does not recognize coalition structure.
Nk=N-2.^(k-1);
vi=v(Nk);
%vi=v(bitset(0,k));
%vnk=v(Nk);
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if sum(vi)<v(N)
   error('sum of lower bound exceeds value of grand coalition! No solution can be found that satisfies the constraints.')
end
%lb=[ra,-Inf]
lb=[];
ub=[vi,Inf];
%ub=[];
if iscell(cs)
   cs=clToMatlab(cs);
else
  cs=double(cs);
end
lcs=length(cs);

S=1:N;
for k=1:n, A1(:,k) = bitget(S,k);end
a1=A1(cs,:);
A1(N+1:N+lcs,:)=-a1;
A1(:,n+1)=1;
A1(N,:)=[];
A1(N:N+lcs-1,end)=0;
A1(cs,end)=0;
A2=sparse(A1);
vc=v;
vc(N)=[];
Bv=v(cs)';
B1=[vc';-Bv];
C=[zeros(1,n),-1];

options.Display='off';
options.Simplex='on';
options.LargeScale='on';
mth1=verLessThan('matlab','24.1.0');
if mth1==0,
    options.Algorithm='dual-simplex-highs';
else
    options.Algorithm='dual-simplex';
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

while 1
  try
    [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
  catch %% old api (before R2022a) with initial value.
   [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
  end
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=find(lambda.ineqlin'>tol);
  bS1(end)=[];
  bA=find(A1(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  it=0:-1:1-n;
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>-tol);
  if exitflag ~= 1
     warning('Prn:Exit','Probably no anti nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     warning('on','all');
     x1=full(x1);
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
end

function [x1, fmin]=weightedCsNucl(clv,cs,pS,tol)
% WEIGHTEDCSNUCL computes a weighted nucleolus of game v w.r.t. coalition structure cs 
% using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=clv.weightedCsNucl(cs,pS,tol)
% Define variables:
%  output:
%  x1        -- A weighted nucleolus of game v w.r.t. cs.
%               (default per capita nucleous).
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  pS       -- A vector of weights of length 2^n-1. (optional)
%  tol      -- Tolerance value. Its default value is set to 10^8*eps. (optional)


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/24/2017        0.9             hme
%   05/25/2024        1.9.2           hme
%                


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
vi=clv.tuvi;

if nargin<3
 pS='';   
 tol=10^6*eps; % Change this value if the solution is not correct.
elseif nargin<4
 tol=10^6*eps;   
end

if iscell(cs)
   cs=clToMatlab(cs);
else
  cs=double(cs);
end
lcs=length(cs);


S=1:N;
for k=1:n, A1(:,k) = -bitget(S,k);end
if isempty(pS)
   mat=-A1';
   clS=ones(1,n)*mat;
   pS=1./clS;
   pS(N)=1;
end
A1=diag(pS)*A1;
a1=A1(cs,:);
A1(N+1:N+lcs,:)=-a1;
A1(:,n+1)=-1;
A1(N,:)=[];
A1(N:N+lcs-1,end)=0;
A1(cs,end)=0;
A2=sparse(A1);
pv=pS.*v;
vc=v;
vc(N)=[];
Bv=vc(cs)';
B1=[-pv';Bv];


C=[zeros(1,n),1];

ra = clv.reasonable_outcome();
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

options.Simplex='on';
options.LargeScale='on';
options.Display='off';
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
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],options);
  catch  
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],options);
  end
  x=xmin;
  x1=x';
  x1(end)=[];
  lambda.ineqlin;
  bS1=find(lambda.ineqlin'>tol);
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
     warning('Prn:Exit','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     x1=full(x1);
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
end

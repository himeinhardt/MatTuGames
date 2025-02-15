function [x1, fmin]=weightedPreNucl_lpp(clv,pS,tol)
% WEIGHTED_PRENUCL_LLP computes a weighted pre-nucleolus of game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=weightedPreNucl_llp(v,pS,tol)
% Define variables:
%  output:
%  x1        -- A weighted pre-nucleolus of game v.
%               (default per capita pre-nucleous).
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  pS       -- A vector of weights of length 2^n-1.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/30/2015        0.6             hme
%   03/29/2015        0.7             hme
%   06/22/2023        1.9.1           hme
%   05/25/2024        1.9.2           hme
%                



if nargin<2
 pS='';   
 tol=10^6*eps; % Change this value if the solution is not correct.
elseif nargin<3
 tol=10^6*eps;   
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;
vi=clv.tuvi;
essQ=clv.tuessQ;
if essQ==0
   error('Game is not essential!');
end

S=1:N;
for k=1:n, A1(:,k) = -bitget(S,k);end
if isempty(pS)
   mat=-A1';
   clS=ones(1,n)*mat;
   pS=1./clS;
   pS(N)=1;
end
A1=diag(pS)*A1;
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
pv=pS.*v;
B1=[-pv';v(N)];
C=[zeros(1,n),1];
Aeq=[ones(1,n),0];
beq=v(N);

ra = clv.reasonable_outcome();
ub=[ra,inf];
x1=[];
lb=[-inf(1,n),-inf];
%ub=[ra,Inf];

opts.Display='off';
opts.Simplex='on';
%opts.ActiveSet='on';
opts.LargeScale='on';
mth1=verLessThan('matlab','24.1.0');
if mth1==0,
    opts.Algorithm='dual-simplex-highs';
else
    opts.Algorithm='dual-simplex';
end
opts.TolFun=1e-10;
opts.TolX=1e-10;
opts.TolRLPFun=1e-10;
%opts
%% for dual-simplex
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);


bA=find(A1(:,end)==0)';
while 1
  try
    [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
  catch %% old api (before R2022a) with initial value.
    [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
  end
  x=xmin';
  x1=x;
  if isempty(x1) == 1
     warning('Prn:Exit0','Probably no pre-nucleolus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  bS1=find(lambda.ineqlin'>0);
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  it=0:-1:1-n;
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
  y=x1;
end

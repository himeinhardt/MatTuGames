function [x1, fmin, pS]=PropNucl(clv,tol)
% PROPNUCL computes the proportional nucleolus of game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=clv.PropNucl(tol)
% Define variables:
%  output:
%  x1        -- The proportional nucleolus of game v.
%  fmin      -- The minmax excess value.
%  pS        -- A weight system.
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
%   05/21/2015        0.7             hme
%                


if nargin<2
 tol=10^6*eps; % Change this value if the solution is not correct.
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
vi=clv.tuvi;

if N==3
  x1=clv.StandardSolution();
  return
end
S=1:N;
snQ=all(sign(v)==-1);
if snQ==0
   pS=1./v;
else
   pS=-1./v;
end
pS(N)=1;
z1=any(isinf(pS));
z2=any(isnan(pS));
    if z1==1 || z2==1
        error('At least one weight is zero!');
    end
for k=1:n, A1(:,k) = -pS.*bitget(S,k);end
A1(N+1,:)=-A1(N,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=A1;
pv=-ones(1,N);
pv(N)=v(N);
B1=[-pv';pv(N)];
C=[zeros(1,n),1];

k=1:n;
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


opts.Display='off';
opts.Simplex='on';
opts.LargeScale='on';
opts.Algorithm='dual-simplex';
opts.TolFun=1e-10;
opts.TolX=1e-10;
opts.TolRLPFun=1e-10;
%% for dual-simplex
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);
bA=find(A1(:,end)==0)';
y=-inf(1,n);
while 1
  [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
  x=xmin';
  x1=x;
  if isempty(x1) == 1
     warning('Prn:Exit0','Probably no proportional nucleolus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  bS1=find(lambda.ineqlin'>tol);
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
     warning('Prn:Exit1','Probably no proportional nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  bA=[bA,bS2];
  y=x1;
end

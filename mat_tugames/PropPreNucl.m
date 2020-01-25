function [x1, fmin, pS]=PropPreNucl(v,tol)
% PROPPRENUCL computes the proportional pre-nucleolus of game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin, pS]=PropPreNucl(v,tol)
% Define variables:
%  output:
%  x1        -- The proportional pre-nucleolus of game v.
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

N=length(v);
[~, n]=log2(N);
if N==3
  x1=StandardSolution(v);
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

ra = reasonable_outcome(v);
ub=[ra,inf];
x1=[];
lb=[-inf(1,n),-inf];

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
     warning('Prn:Exit0','Probably no proportional pre-nucleolus found!')
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
     warning('Prn:Exit1','Probably no proportional pre-nucleolus found!')
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

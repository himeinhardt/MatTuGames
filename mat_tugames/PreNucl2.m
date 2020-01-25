function [x1, fmin]=PreNucl2(v,x1,tol)
% PRENUCL2 computes the pre-nucleolus of game v using the optimization toolbox.
% Use with care, produces round-off errors if simplex method is off.
% Otherwise, very slow.
%
% Usage: [x, fmin]=PreNucl2(v,x1tol)
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of game v.
%  fmin      -- The minmax excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x1       -- starting point for algorithm active-set.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/22/2012        0.3             hme
%   10/21/2014        0.5             hme
%                

N=length(v);
[~, n]=log2(N);

if nargin<3
 tol=10^6*eps; % Change this value if the solution is not correct.
 x1=v(N)*ones(1,n)/n;
elseif nargin<2
 tol=10^6*eps;
end
if N==3
  x1=StandardSolution(v);
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
Aeq=[ones(1,n),0];
beq=v(N);

ra = reasonable_outcome(v);
ub=[ra,Inf];
lb=[-inf(1,n),-Inf];

opts.Display='off';
%opts.Diagnostics='off';
opts.Simplex='on';
%opts.ActiveSet='on';
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
  [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,x1,opts);
  x=xmin;
  x1=x';
  if isempty(x1) == 1
     warning('on','all');
     warning('Prn:Exit0','Probably no pre-nucleolus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  bS1=find(lambda.ineqlin'>tol);
  bA=find(A1(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     warning('on','all');
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
     warning('on','all');
     warning('Prn:Exit1','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     warning('on','all');
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  y=x1;
end

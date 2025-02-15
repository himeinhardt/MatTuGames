function [x1, fmin]=weightedAnti_PreNucl(v,pS,tol)
% WEIGHTEDANTI_PRENUCL computes a weighted anti pre-nucleolus of game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=weightedAnti_PreNucl(v,pS,tol)
% Define variables:
%  output:
%  x1        -- A weighted anti pre-nucleolus of game v.
%               (default per capita pre-nucleous).
%  fmin      -- The minmax excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  pS       -- A vector of weights of length 2^n-1.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/09/2017        0.9             hme
%   06/22/2023        1.9.1           hme
%   05/01/2024        1.9.2           hme
%                

if nargin<2
 pS='';   
 tol=10^6*eps; % Change this value if the solution is not correct.
elseif nargin<3
 tol=10^6*eps;   
end

N=length(v);
[~, n]=log2(N);
S=1:N;
for k=1:n, A1(:,k) = bitget(S,k);end
if isempty(pS)
   mat=-A1';
   clS=ones(1,n)*mat;
   pS=1./clS;
   pS(N)=1;
end
A1=diag(pS)*A1;
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=1;
A1(N:N+1,end)=0;
A2=sparse(A1);
pv=pS.*v;
B1=[pv';-v(N)];
C=[zeros(1,n),-1];
Aeq=[ones(1,n),0];
beq=v(N);

ra = reasonable_outcome(v);
ub=[ra,inf];
x1=[];
lb=[-inf(1,n),-inf];
%ub=[ra,Inf];

opts.Simplex='on';
opts.Display='off';
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
y=-inf(1,n);

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
     warning('Prn:Exit1','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     x1=full(x1); 
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  y=x1;
end

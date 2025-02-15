function [x1,fmin]=PreNucl_mod3(v,tol)
% PRENUCL_MOD3 computes the pre-nucleolus of game v using the optimization toolbox.
% Uses Dual-Simplex (Matlab R2015a).
% 
%
% Usage: [x, fmin]=PreNucl_mod3(v,tol)
%
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of game v.
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
%   08/04/2023        1.9.2           hme
%   04/22/2024        1.9.2           hme
%                


warning('off','all');
if nargin<2
 tol=10^6*eps;
end
%tol=-tol;
mmx(40);
N=length(v);
[~, n]=log2(N);
if N==3
  x1=StandardSolution(v);
  fmin=-inf;
  return
end
% solver parameter
x0=[];
ub=[];
%lb=[];
%ub=[UpperPayoff(v)';inf]; %% produces wrong results.
lb=[smallest_amount(v)';-inf];
% solver parameter
mtv=verLessThan('matlab','9.1.3');
opts.Display='off';
%opts.Diagnostics='off';
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
%% for dual-simplex
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);


S=1:N;
N1=N+1;
A1=zeros(N,n);
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N1,end)=0;
A2=sparse(A1);
v1=v+tol;
B1=[-v1';v(N)];
C=[zeros(n,1);1];
%bA=find(A1(:,end)==0)';
bA=[N,N1];
it=0:-1:1-n;
sS1=[];
eS1=[];
cl=[S,N1];
options.RECT = true;
options.TRANSA = false;
y=-inf(1,n);
%t0=0;
%rk=0;
while 1
%tic;
  try
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
  catch %% old api (before R2022a) with initial value.
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
  end
%  exitflag
%tm=toc;
%t0=tm+t0
  x=xmin;
  x1=x';
  if isempty(x1) == 1
     warning('Prn:Exit0','Probably no pre-nucleolus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  y1=lambda.ineqlin;
  ly1=length(y1);
  if isempty(sS1) && ly1==N1 
     bS1=(find(y1>tol))';
  elseif isempty(sS1) && ly1<N1
     y=zeros(N1,1);
     y(cl)=y1;
     bS1=(find(y>tol))';
  else
     y=zeros(N1,1);
     ccl=cl(ismembc(cl,eS1)==0);
     y(ccl)=y1;
     bS1=(find(y>tol))';
  end
  bS2=bS1(ismembc(bS1,bA)==0);
  if isempty(bS2)==1
     warning('on','all');
     break;
  else
     mbS2=rem(floor(bS2(:)*pow2(it)),2)';
  end
  mS2=rem(floor(bA(:)*pow2(it)),2);
  %mS2(1:2,:)=[];
  tmS2=mS2';
  rk=rank(mS2);
  if exitflag ~=1
     warning('on','all');
     warning('Prn:Exit','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n 
     warning('on','all');
     break;
  else
    if isempty(tmS2)
       CbS2=[];
    elseif rank(tmS2) < rank([tmS2,mbS2]) % no solution
       CbS2=[];	    
    else
       CbS2=linsolve(tmS2,mbS2,options);
    end
    if isempty(CbS2)==0
       mCb3=mmx('mult',tmS2,CbS2);
       cbS2=abs(mCb3-mbS2)<tol;
       cbwz=all(cbS2==1);
       bS2=bS2(cbwz);
    end
  end
%% Checking if coalitions sS are in the span.
%% This has a negative impact on the performance.
%  mC=tmS2\C1;
  sS=S(ismembc(S,bA)==0);
  if rk>2
     C1=rem(floor(sS(:)*pow2(it)),2)';
%  dA = decomposition(tmS2);
     if rank(tmS2) == rank([tmS2,C1])
        mC=linsolve(tmS2,C1,options);
        C2=mmx('mult',tmS2,mC);
        lC=abs(C2-C1)<tol;      
     else % no solution
        mC=[];
     end
  else
     mC=[];
  end
  if isempty(mC)
      cwz=[];
  else
      cwz=all(lC==1);
  end
% Lower this bound below if you do not have 
% enough memory.  
  if n>21
     clear mC C1 mS2 CbS2 tmS2;
  end
  sS1=sS(cwz);
  eS1=unique([sS1,eS1]);
  rS1=ismembc(cl,sS1);
  rbS2=ismembc(cl,bS2);
  A1(rbS2,end)=0;
  B1(rbS2)=B1(rbS2)+fmin;
  A1(rS1,:)=[];
  B1(rS1)=[];
  A2=sparse(A1);
  if isempty(sS1)
  else
     cl=cl(ismembc(cl,sS1)==0);
  end
  bA=unique([bA,bS2,sS1]);
  y=x1;
end

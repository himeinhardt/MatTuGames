function [x1,fmin]=PreNucl_mod4(clv,tol)
% PRENUCL_MOD3 computes the pre-nucleolus of game v using the optimization toolbox.
% Uses Dual-Simplex (Matlab R2015a).
% 
%
% Usage: [x, fmin]=clv.PreNucl_mod3(tol)
%
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of game v.
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/30/2023        1.9.2           hme
%   04/22/2024        1.9.2           hme
%                


warning('off','all');
if nargin<2
 tol=10^6*eps;
end
%tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
if N==3
  x1=clv.StandardSolution();
  fmin=-inf;
  return
end
% solver parameter
x0=[];
ub=[];
%lb=[];
%ub=[clv.UpperPayoff()';inf]; %% produces wrong results.
lb=[clv.smallest_amount()';-inf];
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
eS1=[];
while 1
  try
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
  catch %% old api (before R2022a) with initial value.
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
  end	
  x=xmin;
  x1=x';
  x1(end)=[];
  y1=lambda.ineqlin;
  bS1=(find(y1>tol))';
  bS2=bS1(ismembc(bS1,bA)==0);
  if isempty(bS2)==1
     warning('on','all');
     break;
  else
     mbS2=rem(floor(bS2(:)*pow2(it)),2);
  end
  mS2=rem(floor(bA(:)*pow2(it)),2);
  mS2(1:2,:)=[];
  tmS2=mS2';
  rk=rank(mS2);
  if rk==n
  %   x=(-mS2\B1(bA))';
     break;
  else
    if isempty(mS2)
       CbS2=[];
    elseif rank(tmS2) < rank([tmS2,mbS2']) % no solution
       CbS2=[];       
    else 
       CbS2=mbS2/mS2;
    end
    if isempty(CbS2)==0
       cbS2=(abs(CbS2*mS2-mbS2)<tol)';
       cbwz=all(cbS2==1);
       bS2=bS2(cbwz);
    end
  end
%% Checking if coalitions sS are in the span.
%% This has a negative impact on the performance.
  sS=S(ismembc(S,bA)==0);
  if rk>2
     C1=(rem(floor(sS(:)*pow2(it)),2));
     if rank(tmS2) == rank([tmS2,C1'])     
         mC=C1/mS2;	  
         lC=(abs(mC*mS2-C1)<tol)';
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
  sS1=sS(cwz);
  eS1=[eS1,sS1];
  bS2=[bS2,sS1];
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  bA=unique([bA,bS2]);
end

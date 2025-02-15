function [x1,fmin]=msk_prenucl_mod4(v,tol)
% MSK_PRENUCL_MOD4 computes the pre-nucleolus of game v using mosekmex.
%
% MSK-SOLVER: http://www.mosek.com/
% 
% Usage: [x, alp]=msk_prenucl_mod4(v,tol)
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
%   08/05/2023        1.9.2           hme
%   05/26/2024        1.9.2           hme
%                


warning('off','all');
if nargin<2
 tol=10^6*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);
if N==3
  x1=StandardSolution(v);
  fmin=-inf;
  return
end
% upper bound increases elapsed computation time.
%ra = reasonable_outcome(v);
%ub=[ra,inf]';
% ub=[UpperPayoff(v)';inf]; produces wrong results
lb=[smallest_amount(v)';-inf];
x0=[];
ub=[];
%lb=[];
S=1:N;
N1=N+1;
A1=zeros(N,n);
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
v1=v+tol;
B1=[-v1';v(N)];
C=[zeros(n,1);1];
bA=[N,N1];
it=0:-1:1-n;
eS1=[];
% Define Problem
prob.buc=B1;
prob.a=A2;
c=[zeros(n,1);1];
prob.c=C;
prob.blc=-inf(N+1,1);
prob.blx=lb;
prob.bux=ub;

% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
%param.MSK_IPAR_INTPNT_BASIS   = sc.MSK_OFF;
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex.
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
%param.MSK_DPAR_BASIS_TOL_X = 1.0e-9;
%param.MSK_DPAR_BASIS_TOL_S = 1.0e-9;
%param=[];
eS1=[];
while 1
  [rcode,res] = mosekopt('minimize echo(0)',prob,param);
  sol=res.sol;
  x=sol.bas.xx';
  x1=x;
  x1(end)=[];
  fmin=sol.bas.pobjval;  
  y1=sol.bas.y;  
  bS1=(find(y1<tol))';
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
       cbS2=abs(CbS2*mS2-mbS2 > tol)';
       cbwz=all(cbS2==1);
       bS2=bS2(~cbwz);
    end
  end
%% Checking if coalitions sS are in the span.
%% This has a negative impact on the performance.
  sS=S(ismembc(S,bA)==0);
  if rk>2
     C1=(rem(floor(sS(:)*pow2(it)),2));
     if rank(tmS2) == rank([tmS2,C1'])     
         mC=C1/mS2;       
         lC=abs(mC*mS2-C1 > tol)';
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
     clear mC C1 mS2 CbS2;
  end
  sS1=sS(cwz);
  eS1=[eS1,sS1];
  bS2=[bS2,sS1];
  A1(bS2,end)=0;
  A2=sparse(A1);
  prob.a=A2;  
  B1(bS2)=B1(bS2)+fmin;
  prob.buc=sparse(B1);  
  bA=unique([bA,bS2]);
end

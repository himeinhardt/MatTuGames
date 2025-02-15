function [x1, alp]=msk_cs_nucl(clv,cs,tol)
% MSK_NUCL computes the nucleolus of game v w.r.t. coalition structure cs 
% using mosekmex.
% 
% MSK-SOLVER: http://www.mosek.com/
%
% Usage: [x, alp]=clv.msk_cs_nucl(cs,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game v.
%  alp       -- The minmax excess value.
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
%   09/23/2017        0.9             hme
%                



if nargin<3
 tol=10^8*eps;
end
tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
vi=clv.tuvi;
if N==3
  x1=clv.StandardSolution();
  return
end

S=1:N;

% solver parameter
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
if iscell(cs)
   cs=clToMatlab(cs);
else
  cs=double(cs);
end
lcs=length(cs);

for k=1:n, A1(:,k) = -bitget(S,k);end
a1=A1(cs,:);
A1(N+1:N+lcs,:)=-a1;
A1(:,n+1)=-1;
A1(N,:)=[];
A1(N:N+lcs-1,end)=0;
A2=sparse(A1);
vc=v;
vc(N)=[];
Bv=v(cs)';
B1=[-vc';Bv];

prob.buc=sparse(B1);
c=[zeros(1,n),1];
prob.c=c';
prob.blc=-inf(N+lcs-1,1);
prob.blx=lb';
prob.bux=ub';
% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
%param.MSK_IPAR_INTPNT_BASIS   = 'MSK_ON';
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex (MSK 7).
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
%param.MSK_DPAR_BASIS_TOL_X = 1.0e-9;
%param.MSK_DPAR_BASIS_TOL_S = 1.0e-9;
%param=[];


while 1
  A2=sparse(A1);
  prob.a=A2;
  [rcode,res] = mosekopt('minimize echo(0)',prob,param); 
  sol=res.sol;
  x=sol.bas.xx';
  x1=x;
  x1(end)=[];
  alp=sol.bas.pobjval;
  bS1=(find(sol.bas.y<tol));
%  bS1(end)=[];
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
  posQ=all(wgh>tol);
  if strcmp(sol.bas.solsta,'OPTIMAL') ~= 1
     warning('Prn:ExitB','Probably no nucleolus found!');
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  B1(bS2)=B1(bS2)+alp;
  prob.buc=sparse(B1);
end

function [x1, alp]=msk_nucl(clv,tol)
% MSK_NUCL computes the nucleolus of game v using mosekmex.
% 
% MSK-SOLVER: http://www.mosek.com/
%
% Usage: [x, alp]=msk_nucl(clv,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game v.
%  alp       -- The minmax excess value.
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
%   01/11/2013        0.3             hme
%                



if nargin<2
 tol=10^8*eps;
end
tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
vi=clv.tuvi;
if essQ==0
   error('Game is not essential!')
end
if N==3
  x1=clv.StandardSolution();
  alp=-inf;
  return
end

S=1:N;

% solver parameter
ra = clv.reasonable_outcome();
k=1:n;
vi=v(bitset(0,k));
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
lb=[vi,-Inf];
ub=[ra,Inf];

A1=zeros(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat(1:N,:);
A1(N+1,:)=PlyMat(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
B1=[-v';v(N)];
prob.buc=sparse(B1);
c=[zeros(1,n),1];
prob.c=c';
prob.blc=-inf(N+1,1);
prob.blx=lb';
prob.bux=ub';
% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
param.MSK_IPAR_INTPNT_BASIS   = 'MSK_ON';
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex (MSK 7).
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
%param.MSK_DPAR_BASIS_TOL_X = 1.0e-9;
%param.MSK_DPAR_BASIS_TOL_S = 1.0e-9;
%param=[];
it=0:-1:1-n;

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

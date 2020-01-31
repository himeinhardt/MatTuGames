function [x1, alp]=gurobi_prenucl_llp(v,tol)
% GUROBI_PRENUCL_LLP computes the pre-nucleolus of game v using gurobimex.
% 
% GUROPI-SOLVER: http://www.gurobi.com
%
% Usage: [x, alp]=gurobi_prenucl_llp(v,tol)
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of game v.
%  alp       -- The minmax excess value.
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
%   12/21/2014        0.6             hme
%                



if nargin<2
 tol=10^6*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);
if N==3
  x1=StandardSolution(v);
  return
end
S=1:N;

ra = reasonable_outcome(v);
ub=[ra,Inf];
lb=[-inf(1,n),-Inf];

for k=1:n, A1(:,k) = -bitget(S,k);end
A1(:,end+1)=-1;
A1(N,end)=0;
B1=-v';
objective=[zeros(1,n),1];
s0='<';
ctype=repmat(s0,1,N);
ctype(N)='=';
params.outputflag = 0;
% params.method= 0; % Use primal simplex method.
% params.method= 1; % Use dual simplex method.
% params.method= 2; % Use barrier method.
params.method= 3; % Use concurrent.
% params.method= 4; % Use deterministic concurrent.
params.FeasibilityTol = 1e-9;
params.OptimalityTol = 1e-9;
params.BarConvTol = 1e-10;
%params.TimeLimit = 3000;
it=0:-1:1-n;
bA=find(A1(:,end)==0);
while 1
  A2=sparse(A1);
  model.A=A2;
  model.obj=objective;
  model.rhs = B1;
  model.sense = ctype;
  model.vtype = 'C';
  model.modelsense = 'min';
  model.lb=lb';
  model.ub=ub';
  result = gurobi(model,params);
  if strcmp(result.status,'INFEASIBLE')
     warning('Prn:ExitA','Probably no pre-nucleolus found!');
     break;
  end
  x=result.x;
  x1=x';
  x1(end)=[];
  alp=result.objval;
  bS1=(find(result.pi<tol));
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  bA=[bA;bS2];
  mS2=rem(floor(bA(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+alp;
  if rk==n 
     x=(-mS2\B1(bA))';
     break;
  end
  A1(bS2,end)=0;
  ctype(bS2)='=';
end

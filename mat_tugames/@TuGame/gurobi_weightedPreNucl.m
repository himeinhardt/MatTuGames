function [x1, alp]=gurobi_weightedPreNucl(clv,pS,tol)
% GUROBI_WEIGHTEDPRENUCL computes a weighted pre-nucleolus of game v using gurobimex.
% 
% GUROPI-SOLVER: http://www.gurobi.com
%
% Usage: [x, alp]=gurobi_weightedPreNucl(clv,pS,tol)
% Define variables:
%  output:
%  x1        -- A weighted pre-nucleolus of game v.
%               (default per capita pre-nucleous).
%  alp       -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  pS       -- A vector of weights of length 2^n-1.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/03/2014        0.5             hme
%                


if nargin<2
 pS='';   
 tol=10^6*eps; % Change this value if the solution is not correct.
elseif nargin<3
 tol=10^6*eps;   
end
tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;


ra = clv.reasonable_outcome();
ub=[ra,Inf];
lb=[-inf(1,n),-Inf];

S=1:N;
for k=1:n, A1(:,k) = -bitget(S,k);end
if isempty(pS)
   mat=-A1';
   clS=ones(1,n)*mat;
   pS=1./clS;
   pS(N)=1;
end
A1=diag(pS)*A1;
A1(:,end+1)=-1;
A1(N,end)=0;
pv=pS.*v;
B1=-pv';
objective=[zeros(1,n),1];
s0='<';
ctype=repmat(s0,1,N);
ctype(N)='=';
% params.method= 0; % Use primal simplex method.
% params.method= 1; % Use dual simplex method.
% params.method= 2; % Use barrier method.
params.method= 3; % Use concurrent.
% params.method= 4; % Use deterministic concurrent.
params.FeasibilityTol = 1e-9;
params.OptimalityTol = 1e-9;
params.BarConvTol = 1e-10;
params.TimeLimit = 1000;


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
  params.outputflag = 0;
  params.method= 1; % Use dual simplex method.
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
  if strcmp(result.status,'OPTIMAL') ~= 1
     warning('Prn:ExitB','Probably no pre-nucleolus found!');
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  ctype(bS2)='=';
  B1(bS2)=B1(bS2)+alp;
end

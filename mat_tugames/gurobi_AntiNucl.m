function [x1, alp]=gurobi_AntiNucl(v,tol)
% GUROBI_ANTINUCL computes the anti nucleolus of game v using gurobimex.
% 
% GUROPI-SOLVER: http://www.gurobi.com
%
% Usage: [x, alp]=gurobi_AntiNucl(v,tol)
% Define variables:
%  output:
%  x1        -- The anti nucleolus of game v.
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
%   02/09/2017        0.9             hme
%                


warning('off','all');
if nargin<2
 tol=10^6*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);

S=1:N;

% solver parameter
ra = smallest_amount(v);
k=1:n;
Nk=N-2.^(k-1);
vi=v(Nk);
%vi=v(bitset(0,k));
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if sum(vi)<v(N)
   error('sum of lower bound exceeds value of grand coalition! No solution can be found that satisfies the constraints.')
end
lb=[-ra,-Inf];
ub=[vi,Inf];

for k=1:n, A1(:,k) = bitget(S,k);end
A1(:,end+1)=1;
A1(N,end)=0;
B1=v';
objective=[zeros(1,n),-1];
s0='<';
ctype=repmat(s0,1,N);
ctype(N)='=';
params.outputflag = 0;
% params.method= 0; % Use primal simplex method.
% params.method= 1; % Use dual simplex method.
params.method= 2; % Use barrier method.
% params.method= 3; % Use concurrent.
% params.method= 4; % Use deterministic concurrent.
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
  result = gurobi(model,params);
  if strcmp(result.status,'INFEASIBLE')
     warning('on','all');
     warning('Prn:ExitA','Probably no anti nucleolus found!');
     break;
  elseif strcmp(result.status,'INF_OR_UNBD')
     warning('on','all');
     warning('Prn:ExitA','Probably no anti nucleolus found!');
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
     warning('on','all');
     warning('Prn:ExitB','Probably no anti nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  ctype(bS2)='=';
  B1(bS2)=B1(bS2)+alp;
end

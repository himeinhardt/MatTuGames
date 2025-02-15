function [x1, alp]=gurobi_cfr_nucl(v,F,tol)
% GUROBI_CFR_NUCL computes the nucleolus of game v with coalition formation restrictions 
% using gurobimex.
% 
%
% Source: Granot et al. (1978), Characterization sets for the nucleolus. IJGT.
%
%
% GUROPI-SOLVER: http://www.gurobi.com
%
% Usage: [x, alp]=gurobi_cfr_nucl(v,F)
% Define variables:
%  output:
%  x1        -- The nucleolus of game vF.
%  alp       -- The minmax excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  F        -- For instance, a characterization set for the nucleolus.
%              F must contain the grand coalition N. 
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/15/2017        0.9             hme
%                



if nargin<2
   error('A collection of sets F is required!');
elseif nargin<3
 tol=10^6*eps;
end

tol=-tol;

N=length(v);
[~, n]=log2(N);
k=1:n;
si=bitset(0,k);
%% Addining the singleton coalitions to F, 
%% since v is definded over F ans si.
F=unique([F,si]);
lf=length(F);


% solver parameter
ra = reasonable_outcome(v);
k=1:n;
vi=v(si);
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

%% F should contain the grand coalition for defining vF.
S=1:N;
lfNq=F(end)~=N;
if lfNq
   F(end+1)=N;
   lf=lf+1;
end
CS=S(ismember(S,F)==0);
vF=v;
vF(CS)=[];
if lfNq
   vF(end)=0;
end


for k=1:n, A1(:,k) = -bitget(F,k);end
A1(:,end+1)=-1;
A1(lf,end)=0;
B1=-vF';
objective=[zeros(1,n),1];
s0='<';
ctype=repmat(s0,1,lf);
ctype(lf)='=';
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
     warning('Prn:ExitA','Probably no nucleolus found!');
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
  bA(end)=[];
  bA=F(bA);
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>tol);
  if strcmp(result.status,'OPTIMAL') ~= 1
     warning('Prn:ExitB','Probably no nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  ctype(bS2)='=';
  B1(bS2)=B1(bS2)+alp;
end

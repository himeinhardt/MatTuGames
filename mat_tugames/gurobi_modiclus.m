function [x1, alp]=gurobi_modiclus(v,tol)
% GUROBI_MODICLUS computes the modiclus of game v using gurobimex.
% 
% GUROPI-SOLVER: http://www.gurobi.com
%
% Usage: [x, alp]=gurobi_modiclus(v,tol)
% Define variables:
%  output:
%  x1        -- The modiclus of game v.
%  alp       -- The minmax bi-excess value.
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
%   12/19/2019        0.9             hme
%                



if nargin<2
 tol=10^6*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);

S=0:N;
N1=N+1;
n1=2*n;
N2=(N1)^2-N1;
for k=1:n, A1(:,k) = -bitget(S,k);end
ve=[0,v];
B1=zeros(N2,1);
A2=zeros(N2,n);
cl=zeros(1,N2);
ii=1;
for k=1:N1
    for jj =1:N1
        if k ~= jj
           if k>1 && jj >1 
              cl(ii)=(k-1)+(jj-1)*N1;
           elseif k==1 && jj >1
              cl(ii)=N1*(jj-1);
           elseif k>1 && jj==1
              cl(ii)=k-1;
           end
           A2(ii,:) = A1(k,:)-A1(jj,:);       
           B1(ii) = ve(jj)-ve(k);         
           ii = ii+1;
        end
    end
end


ra = reasonable_outcome(v);
ub=[ra,Inf];
%ub=inf(1,n+1);
lb=[-inf(1,n),-Inf];
cl(ii)=2^n1-1;
A2(ii,:)=-ones(1,n);
A2(:,end+1)=-1;
A2(ii,end)=0;
%B1=[-vs;v(N)+dv(N)];
B1=[B1;-v(N)];
objective=[zeros(1,n),1];
s0='<';
ctype=repmat(s0,1,N2+1);
ctype(N2+1)='=';
params.outputflag = 0;
% params.method= 0; % Use primal simplex method.
% params.method= 1; % Use dual simplex method.
% params.method= 2; % Use barrier method.
params.method= 3; % Use concurrent.
% params.method= 4; % Use deterministic concurrent.
params.FeasibilityTol = 1e-9;
params.OptimalityTol = 1e-9;
params.BarConvTol = 1e-10;
% params.Threads=8;
%params.TimeLimit = 3000;
it=0:-1:1-n1;
bA=find(A2(:,end)==0);
while 1
  A2=sparse(A2);
  model.A=A2;
  model.obj=objective;
  model.rhs = B1;
  model.sense = ctype;
  model.vtype = 'C';
  model.modelsense = 'min';
  model.lb=lb';
  model.ub=ub';
  result = gurobi(model,params);
  if strcmp(result.status,'OPTIMAL')==0
     warning('Prn:ExitA','Probably no modiclus found!');
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
  S2=cl(bA);
  mS2=rem(floor(S2(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+alp;
  if rk==n1
     x=(-mS2\B1(bA))';
     break;
  end
  A2(bS2,end)=0;
  ctype(bS2)='=';
end

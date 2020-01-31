function [x1, alp]=msk_modiclus(v,tol)
% MSK_MODICLUs computes the modiclus of game v using mosekmex.
% 
% MSK-SOLVER: http://www.mosek.com/
%
% Usage: [x, alp]=msk_modiclus(v,tol)
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
%   12/19/2017        0.9             hme
%                



if nargin<2
 tol= 10^8*eps; % Change this value if the solution is not correct.
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


% upper bound increases elapsed computation time.
ra = reasonable_outcome(v);
ub=[ra,inf]';
lb=-inf(n+1,1);


A2(N2+1,:)=-ones(1,n);
A2(:,end+1)=-1;
A2(N2+1,end)=0;
B1(N2+1)=-v(N);
%A1(N+1,end)=0;
prob.buc=B1;
c=[zeros(n,1);1];
prob.c=c;
prob.blc=-inf(N2,1);
prob.blc(N2+1)=-v(N);
prob.blx=lb;
%prob.bux=ub;
% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
%param.MSK_IPAR_INTPNT_BASIS   = sc.MSK_OFF;
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex (MSK 7).
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
%param.MSK_DPAR_BASIS_TOL_X = 1.0e-9;
%param.MSK_DPAR_BASIS_TOL_S = 1.0e-9;
%param=[];
it=0:-1:1-n1;

while 1
  A2=sparse(A2);
  prob.a=A2;
  [rcode,res] = mosekopt('minimize echo(0)',prob,param);
  sol=res.sol;
  x=sol.bas.xx';
  x1=x;
  x1(end)=[];
  alp=sol.bas.pobjval;
  bS1=(find(sol.bas.y<tol));
  AA=A2;
  AA(N2+1,:)=[];
  bA=find(AA(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  S2=cl(bA);
  mS2=rem(floor(S2(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+alp;
  if rk==n1 
     x=(-mS2\B1(bA))';
     break;
  end
  A2(bS2,end)=0;
  prob.buc=sparse(B1);
end

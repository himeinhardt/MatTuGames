function [x1, fmin]=Cost_PreNucl(clv,tol)
% COST_PRENUCL computes the pre-nucleolus of a cost game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=Cost_PreNucl(clv,tol)
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of a cost game v.
%  fmin      -- The maxmin excess value.
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
%   03/25/2021        1.9             hme
%   05/25/2024        1.9.2           hme
%                


if nargin<2
 tol=10^6*eps; % Change this value if the solution is not correct.
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayer
if N==3
  x1=clv.StandardSolution();
  return
end
S=1:N;
A1=zeros(N,n);
for k=1:n, A1(:,k) = bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[v';-v(N)];
C=[zeros(1,n),-1];
%Aeq=[ones(1,n),0];
beq=v(N);

%ra = clv.reasonable_outcome();
ub=[inf(1,n),inf];
x1=[];
lb=[-inf(1,n),-inf];
%ub=[ra,Inf];

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
bA=find(A1(:,end)==0)';
y=-inf(1,n);
it=0:-1:1-n;
while 1
  try
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
  catch
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
  end   
  x=xmin';
  x1=x;
  if isempty(x1) == 1
     warning('Prn:Exit0','Probably no pre-nucleolus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  bS1=find(lambda.ineqlin'>tol);
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>-tol);
  if exitflag ~= 1
     warning('Prn:Exit1','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     x1=full(x1); 
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  bA=[bA,bS2];
  y=x1;
end

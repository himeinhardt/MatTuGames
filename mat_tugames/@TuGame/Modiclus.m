function [x1, fmin]=Modiclus(clv,tol)
% MODICLUS computes the modiclus of game v using the optimization toolbox.
% Requires Dual-Simplex (Matlab R2015a).
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: [x, fmin]=clv.Modiclus(tol)
% Define variables:
%  output:
%  x1        -- The modiclus of game v.
%  fmin      -- The minmax bi-excess value.
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/18/2017        0.9             hme
%   04/22/2024        1.9.2           hme
%                


if nargin<2
 tol=10^6*eps; % Change this value if the solution is not correct.
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;


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
Aeq=[ones(1,n),0];
A2(:,end+1)=-1;
Beq=v(N);
C=[zeros(1,n),1];
ra = clv.reasonable_outcome();
ub=[ra,inf];
x1=[];
lb=-inf(1,n+1);
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
y=-inf(1,n);
it=0:-1:1-n1;
bA=find(A2(:,end)==0)';
while 1
  A2=sparse(A2);
  try
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,Aeq,Beq,lb,ub,opts);
  catch %% old api (before R2022a) with initial value.
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,Aeq,Beq,lb,ub,[],opts);
  end
  x=xmin';
  x1=x;
  if isempty(x1) == 1
     warning('Prn:Exit0','Probably no modiclus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  bS1=find(lambda.ineqlin'>tol);
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  bA=[bA,bS2];
  S2=cl(bA);
  mS2=rem(floor(S2(:)*pow2(it)),2); %% For the balanced set, we can probalby do better!!!
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+fmin;
  if rk==n1
     x=(-mS2\B1(bA))';
     break;
  end
  A2(bS2,end)=0;
  y=x1;
end

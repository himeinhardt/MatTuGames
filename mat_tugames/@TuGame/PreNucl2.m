function [x1, fmin]=PreNucl2(clv,x1,tol)
% PRENUCL2 computes the pre-nucleolus of game v using the optimization toolbox.
% Use with care, produces large round-off errors if simplex method is off.
% Otherwise, very slow.
%
% Usage: [x, fmin]=PreNucl2(clv,x1,tol)
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of game clv.
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  x1       -- starting point
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/27/2013        0.3             hme
%   05/26/2024        1.9.2           hme
%                



v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;
vi=clv.tuvi;

if nargin<3
 tol=10^6*eps; % Change this value if the solution is not correct.
 x1=v(N)*ones(1,n)/n;
elseif nargin<2
 tol=10^6*eps;
end



S=1:N;
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[-v';v(N)];
bS2=[];
C=[zeros(1,n),1];

ra = reasonable_outcome(v);
ub=[ra,Inf];

if strcmp(gt,'cv')
 lb=[vi,-Inf];
else
 lb=[];
end
% produces large round-off errors.
%opts=optimset('TolFun',1e-10,'TolX',1e-10,'MaxIter',128);
%opts=optimset('Simplex','off','LargeScale','on','MaxIter',128);

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


while 1
  try
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
  catch %% old api (before R2022a) with initial value.
     [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,x1,opts);
  end
  x=xmin;
  x1=x';
  if isempty(x1) == 1
     warning('Prn:Exit0','Probably no pre-nucleolus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  lambda.ineqlin;
  bS1=find(lambda.ineqlin'>tol);
  bS1(end)=[];
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
  posQ=all(wgh>-tol);
  if exitflag ~= 1
     warning('Prn:Exit1','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  y=x1;
end

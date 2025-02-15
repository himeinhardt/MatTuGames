function [x1,fmin]=cplex_AntiPreNucl_mod3(clv,tol)
% CPLEX_ANTIPRENUCL_MOD3 computes the anti pre-nucleolus of game v using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
% 
%
% Usage: [x, alp]=clv.cplex_AntiPreNucl_mod3(tol)
% Define variables:
%  output:
%  x1        -- The anti pre-nucleolus of game v.
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
%   04/06/2020        1.9             hme
%                


warning('off','all');
if nargin<2
 tol=10^6*eps;
end
%tol=-tol;
mmx(32);
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
% solver parameter
x0=[];
%ub=[];
%lb=[];
ub=[UpperPayoff(v)';inf];
lb=[smallest_amount(v)';-inf];
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else 
%  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
%  options.MaxIter='128';
  options.display='off';
  options.largescale='on';
  options.algorithm='dual-simplex';
  options.tolfun=1e-10;
  options.tolx=1e-10;
  options.tolrlpfun=1e-10;
  %%%% for dual-simplex
  % opts.MaxTime=9000;
  options.preprocess='none';
  options.tolcon=1e-8;
  options.maxiter=10*(N+n);
  options.threads=3;
end

S=1:N;
N1=N+1;
A1=zeros(N,n);
for k=1:n, A1(:,k) = bitget(S,k);end
A1(N1,:)=-A1(end,:);
A1(:,end+1)=1;
A1(N:N1,end)=0;
A2=sparse(A1);
v1=v+tol;
B1=[v1';-v(N)];
C=[zeros(n,1);-1];
%bA=find(A1(:,end)==0)';
bA=[N,N1];
it=0:-1:1-n;
sS1=[];
eS1=[];
cl=[S,N1];
opts.RECT = true;
opts.TRANSA = false;
y=-inf(1,n);
%t0=0;
%rk=0;
while 1
%tic;
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  if isempty(xmin)
     ub=[];
     [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
%  exitflag
%tm=toc;
%t0=tm+t0
  x=xmin;
  x1=x';
  if isempty(x1) == 1
     warning('APrn:Exit0','Probably no anti pre-nucleolus found!')
     x1=y;
     break;
  end
  x1(end)=[];
  y1=lambda.ineqlin;
  ly1=length(y1);
  if isempty(sS1) && ly1==N1 
     bS1=(find(y1>tol))';
  elseif isempty(sS1) && ly1<N1
     y=zeros(N1,1);
     y(cl)=y1;
     bS1=(find(y>tol))';
  else
     y=zeros(N1,1);
     ccl=cl(ismembc(cl,eS1)==0);
     y(ccl)=y1;
     bS1=(find(y>tol))';
  end
  bS2=bS1(ismembc(bS1,bA)==0);
  if isempty(bS2)==1
     warning('on','all');
     break;
  else
     mbS2=rem(floor(bS2(:)*pow2(it)),2)';
  end
  mS2=rem(floor(bA(:)*pow2(it)),2);
  %mS2(1:2,:)=[];
  tmS2=mS2';
  rk=rank(mS2);
  if exitflag ~=1
     warning('on','all');
     warning('APrn:Exit','Probably no anti pre-nucleolus found!')
     break; 
  elseif rk==n 
     warning('on','all');
     break;
  else
    if isempty(tmS2)
       CbS2=[];
    else
       CbS2=linsolve(tmS2,mbS2,opts);
    end
    if isempty(CbS2)==0
       mCb3=mmx('mult',tmS2,CbS2);
       cbS2=abs(mCb3-mbS2)<tol;
       cbwz=all(cbS2==1);
       bS2=bS2(~cbwz);
    end
  end
%% Checking if coalitions sS are in the span.
%% This has a negative impact on the performance.
  sS=S(ismembc(S,bA)==0);
  C1=rem(floor(sS(:)*pow2(it)),2)';
%  dA = decomposition(tmS2);
  mC=linsolve(tmS2,C1,opts);
%  mC=tmS2\C1;
  if rk>2
      C2=mmx('mult',tmS2,mC);
      lC=abs(C2-C1)<tol;
  else
     mC=[];
  end
  if isempty(mC)
      cwz=[];
  else
      cwz=all(lC==1);
  end
  sS1=sS(cwz);
  eS1=unique([sS1,eS1]);
  rS1=ismembc(cl,sS1);
  rbS2=ismembc(cl,bS2);
  A1(rbS2,end)=0;
  B1(rbS2)=B1(rbS2)+fmin;
  A1(rS1,:)=[];
  B1(rS1)=[];
  A2=sparse(A1);
  if isempty(sS1)
  else
     cl=cl(ismembc(cl,sS1)==0);
  end
  bA=unique([bA,bS2,sS1]);
  y=x1;
end

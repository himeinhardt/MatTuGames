function [x1,fmin]=cplex_AntiPreNucl_mod(clv,tol)
% CPLEX_ANTIPRENUCL_MOD computes the anti pre-nucleolus of game v using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
% 
%
% Usage: [x, alp]=clv.cplex_AntiPreNucl_mod(tol)
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

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
% solver parameter
x0=[];
ub=[];
lb=[];
%cplex = Cplex('null');
%str=cplex.getVersion;
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
  options.Preprocess='none';
  options.tolcon=1e-6;
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
%ov=ones(1,n);
%Aeq=ones(1,n+1);
%beq=v(N);
eS1=[];
while 1
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
%% We try this to overcome numerical issues!!
%exitflag
  if exitflag~=1
      ub=[UpperPayoff(v)';inf];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if exitflag~=1
      ub=[UpperPayoff(v)';inf];
      lb=[smallest_amount(v)';-inf];
      B1=[-v1';v(N)+tol];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if exitflag~=1
      %ub=[UpperPayoff(v)';inf];
      %lb=[smallest_amount(v)';-inf];
      v1=v-tol;
      B1=[-v1';v(N)];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if exitflag~=1
      B1=[-v1';v(N)-tol];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if exitflag~=1
      B1=[-v1';v(N)+tol];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if  exitflag~=1
      warning('on','all');
      warning('APrn:Exit','Probably no anti pre-nucleolus found!');
      x1=[];
      break;
  end
  x=xmin;
  x1=x';
  x1(end)=[];
  y1=lambda.ineqlin;
  bS1=(find(y1>tol))';
  bS2=bS1(ismembc(bS1,bA)==0);
  if isempty(bS2)==1
     warning('on','all');
     break;
  else
     mbS2=rem(floor(bS2(:)*pow2(it)),2)';
  end
  mS2=rem(floor(bA(:)*pow2(it)),2);
  mS2(1:2,:)=[];
  tmS2=mS2';
  rk=rank(mS2);
  if rk==n
  %   x=(-mS2\B1(bA))';
     break;
  else
    if isempty(tmS2)
       CbS2=[];
    else 
       CbS2=linsolve(tmS2,mbS2);
    end
    if isempty(CbS2)==0
       cbS2=abs(tmS2*CbS2-mbS2)<tol;
       cbwz=all(cbS2==1);
       bS2=bS2(~cbwz);
    end
  end
%% Checking if coalitions sS are in the span.
%% This has a negative impact on the performance.
  if rk>2
    sS=S(ismembc(S,bA)==0);
    C1=rem(floor(sS(:)*pow2(it)),2)';
%% Needs for n=22 about 7 TB memory.!!!
    mC=tmS2\C1;
    tm=tmS2*mC;
    lC=abs(tm-C1)<tol;
  else
     mC=[];
  end
  if isempty(mC)==0
      cwz=all(lC==1);
      sS1=sS(cwz);
      eS1=[eS1,sS1];
      bS2=[bS2,sS1];
  end  
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  bA=unique([bA,bS2]);
end

function [x1, fmin]=cplex_prenucl(v,tol)
% CPLEX_PRENUCL computes the pre-nucleolus of game v using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
% 
%
% Usage: [x, alp]=cplex_prenucl(v,tol)
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of game v.
%  fmin      -- The minmax excess value.
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
%   12/21/2012        0.3             hme
%   10/15/2013        0.5             hme
%   02/24/2018        0.9             hme
%   03/16/2019        1.0             hme
%   04/04/2020        1.9             hme
%                


warning('off','all');
if nargin<2
 tol=10^6*eps;
end
%tol=-tol;

N=length(v);
[~, n]=log2(N);
if N==3
  x1=StandardSolution(v);
  return
end
% solver parameter
ub=[];
lb=[];
x0=[];
%cplex = Cplex('null');
%str=cplex.getVersion;
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else 
%  options = cplexoptimset('Algorithm','primal','Display','off');
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
  options.display='off';
  options.threads=3;
end

S=1:N;
A1=zeros(N,n);
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
v1=v+tol;
B1=[-v1';v(N)];
C=[zeros(n,1);1];
bA=find(A1(:,end)==0)';
it=0:-1:1-n;
while 1
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
%% We try this to overcome numerical issues!!
%%exitflag
  if exitflag~=1
      B1=[-v1';v(N)-tol];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if exitflag~=1
      B1=[-v1';v(N)+tol];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if exitflag~=1
      ub=[UpperPayoff(v)';inf];
      [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  end
  if exitflag~=1
      %ub=[UpperPayoff(v)';inf];
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
  x=xmin;
  x1=x';
  x1(end)=[];
  y1=lambda.ineqlin;
  bS1=(find(y1>tol))';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     warning('on','all');
     break;
  end
  mS2=rem(floor(bA(:)*pow2(it)),2);
  mS2(1:2,:)=[];
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
%  wgh=pinv(tmS2)*ov';
  wgh=(y1+1);
  posQ=all(wgh>-tol);
  if exitflag ~=1
     warning('on','all');
     warning('Prn:Exit','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     warning('on','all');
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  bA=[bA,bS2];
end

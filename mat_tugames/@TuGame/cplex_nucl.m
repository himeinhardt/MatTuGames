function [x1, fmin]=cplex_nucl(clv,tol)
% CPLEX_NUCL computes the nucleolus of game v using cplexmex.
% 
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.8.0 and higher)
%
%
% Usage: [x, alp]=cplex_nucl(clv,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game v.
%  fmin      -- The minmax excess value.
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
%   12/28/2012        0.3             hme
%   10/15/2013        0.5             hme
%   02/24/2018        0.9             hme
%   04/01/2019        1.0             hme
%                



if nargin<2
 tol=10^8*eps;
end
%tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
vi=clv.tuvi;
if essQ==0
   error('Game is not essential!')
end
if N==3
  x1=clv.StandardSolution();
  return
end

% solver parameter
ra = clv.reasonable_outcome();
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
lb=[vi,-Inf]';
ub=[ra,Inf]';

x0=[];
warning('off','all');
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else 
  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
  options.LargeScale='on';
  options.Algorithm='dual-simplex';
  options.TolFun=1e-10;
  options.TolX=1e-10;
  options.TolRLPFun=1e-10;
  %%%% for dual-simplex
  % opts.MaxTime=9000;
  options.Preprocess='none';
  options.TolCon=1e-6;
  options.MaxIter=10*(N+n);
end
warning('on','all');
S=1:N;
A1=zeros(N,n);
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[-v';v(N)];
C=[zeros(n,1);1];
it=0:-1:1-n;
while 1
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(lambda.ineqlin>tol))';
  bS1(end)=[];
  bA=find(A1(:,end)==0)';
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
  if exitflag ~=1
     warning('Prn:Exit','Probably no nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
end

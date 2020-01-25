function [x1, fmin]=cplex_nucl2(clv,x1,tol)
% CPLEX_NUCL computes the nucleolus of game v from a starting point using glpkmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.5.1 and higher)
% 
%
% Usage: [x, fmin]=cplex_nucl(v,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game v.
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  x1       -- starting point.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/13/2015        0.7             hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
vi=clv.tuvi';

if nargin<2
 x1=v(N)*ones(1,n)/n;
 tol=10^6*eps;
elseif nargin<3
 tol=10^6*eps;
end

if essQ==0
   error('Game is not essential!')
end
if N==3
  x1=clv.StandardSolution();
  return
end


% solver parameter
ra = clv.reasonable_outcome;
k=1:n;
vi=v(bitset(0,k));
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if sum(vi)>v(N)
   error('sum of lower bound exceeds value of grand coalition! No solution can be found that satisfies the constraints.')
end
lb=[vi,-Inf]';
ub=[ra,Inf]';

warning('off','all');
options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
warning('on','all');
S=1:N;
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[-v';v(N)];
C=[zeros(n,1);1];

while 1
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x1,options);
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
  it=0:-1:1-n;
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

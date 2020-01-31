function [fmin,x,A1,B1,exitflag]=cplex_LeastCore(clv,tol)
% CPLEX_LEASTCORE computes the least core of game using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.8.0 and higher)
%
%
% Usage: [epsv, x]=cplex_LeastCore(v,tol)
% Define variables:
%  output:
%  fmin         -- Critical value of the least core.
%  x            -- An allocation that satisfies the least-core constraints.
%  A1           -- least-core matrix.
%  B1           -- boundary vector.
%  exitflag     -- Exit flag of the LP. Returns 1 for a successful termination.
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
%   06/03/2015        0.7             hme
%   02/24/2018        0.9             hme
%                


if nargin<2
 tol=10^8*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat(1:N,:);
A1(N+1,:)=PlyMat(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[-v';v(N)];
C=[zeros(n,1);1];
ub=[ones(n,1)*v(N);Inf];
warning('off','all');
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  opts = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else 
  opts = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
end
warning('on','all');
[xmin,fmin,exitflag,output,lambda]=cplexlp(C,A2,B1,[],[],[],ub,[],opts);
x=xmin';
x(end)=[];

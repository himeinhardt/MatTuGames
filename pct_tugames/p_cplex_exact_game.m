function [v_e xm]=p_cplex_exact_game(v,tol);
% P_CPLEX_EXACT_GAME computes the exact game from game v using cplexmex
% and Matlab's PCT.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.8.0 and higher)
% 
%  Usage: [v_e xm status]=p_cplex_exact_game(v,tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to -10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/25/2016        0.8             hme
%   02/24/2018        0.9             hme
%                
if nargin<2
 tol=10^6*eps;
end

N=length(v);
[~,n]=log2(N);

if CddCoreQ(v,tol)==0
  error('Core is empty!');
 else
end

k=1:n;
vi=v(bitset(0,k));

S=1:N;
parfor k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat;
A1(N+1,:)=PlyMat(end,:);
A2=sparse(A1);
B1=[-v';v(N)];
C=[zeros(n,1);1];
ub=inf(1,n);
lb=vi;
warning('off','all');
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  opts = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else
  opts = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
end
warning('on','all');

parfor S=1:N
    C=PlyMat(S,:)';
    [xmin,fmin,exitflag]=cplexlp(C,A2,B1,[],[],lb,ub,[],opts);
    xm(S,:)=xmin';
    v_e(S)=fmin;
    if exitflag ~= 1 
       warning('Solution not optimal!')
    end
end    



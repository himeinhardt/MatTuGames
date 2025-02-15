function [alpv, v]=cplex_alphaVector(v,pt,tol)
% CPLEX_ALPHAvECTOR computes recursively an alpha vector from the positive core of 
% game v using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
%
% Usage: [alpv, v]=cplex_alphaVector(v,pt,tol)
% Define variables:
%  output:
%  alpv     -- Alpha vector computed from the positive core.
%  v        -- A derived Tu-Game v of length 2^n-1.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  pt       -- A partition of N by their unique integers representation.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/02/2015        0.6             hme
%   02/24/2018        0.9             hme
%   04/04/2020        1.9             hme
%                


if nargin<3
 tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
% Determining the lower bounds
S=1:N;
%k=1:n;
%vi=v(bitset(0,k));
%lb=vi';
lb=[];

% Determining the pre-nucleolus
try
   x=cplex_prenucl(v);
catch
   x=PreNucl(v);
end
exv=max(excess(v,x),0);
ad=additive_game(x);
cex=S(exv>tol);
v(cex)=ad(cex);

% Constructing the LP
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(N,:);
A2=sparse(A1);
lpt=length(pt);
alpv=zeros(1,lpt);
cmat=-A1(pt,:);
ub=[];
x0=[];

warning('off','all');
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else
%  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
  options.largescale='on';
  options.algorithm='dual-simplex';
  options.tolfun=1e-10;
  options.tolx=1e-10;
  options.tolrlpfun=1e-10;
  %%%% for dual-simplex
  % opts.MaxTime=9000;
  options.preprocess='none';
  options.tolcon=1e-6;
  options.maxiter=128;
  options.display='off';
  options.threads=3;
end
warning('on','all');
%options

for ii=1:lpt
    B1=[-v';v(N)];
    C=cmat(ii,:);
    [x,fmin,exitflag,output,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
    if exitflag<0
       exitflag
       error('Optimization terminated not successfully')
    end
    alpv(ii)=cmat(ii,:)*x;
    v(pt(ii))=alpv(ii);
end

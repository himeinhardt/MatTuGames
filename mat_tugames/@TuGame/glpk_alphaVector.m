function [alpv, v]=glpk_alphaVector(clv,pt,tol)
% GLPK_ALPHAvECTOR computes recursively an alpha vector from the positive core of game v 
% using glpkmex.
%
% http://www.gnu.org/software/glpk/glpk.html
%
% Usage: [alpv, v]=clv.glpk_alphaVector(pt,tol)
% Define variables:
%  output:
%  alpv     -- Alpha vector computed from the positive core.
%  v        -- A derived Tu-Game v of length 2^n-1.
%
%  input:
%  clv      -- TuGame class object.
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
%                


if nargin<3
 tol=10^6*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
% Determining the lower bounds. 
S=1:N;
%k=1:n;
%vi=v(bitset(0,k));
%lb=vi';
lb=[];
% Determining the pre-nucleolus
try
   x=clv.glpk_prenucl();
catch
   x=clv.PreNucl;
end
exv=max(clv.excess(x),0);
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

ctype=[];
vartype=[];
s=1; % minimization problem 
param.lpsolver=1; % simplex method

for ii=1:lpt
    B1=[-v';v(N)];
    C=cmat(ii,:);
    [xmin, fmin, status, extra] = glpk(C,A2,B1,lb,ub,ctype,vartype,s,param);
    if status ~=5
       status
       error('Optimization terminated not successfully')
    end
    alpv(ii)=fmin;
    v(pt(ii))=alpv(ii);
end

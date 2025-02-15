function [alpv, v]=msk_alphaVector(v,pt,tol)
% MSK_ALPHAvECTOR computes recursively an alpha vector from the positive core of game v 
% using mosekmex..
%
%  MSK-SOLVER: http://www.mosek.com/
%
%
% Usage: [alpv, v]=msk_alphaVector(v,pt,tol)
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
%                


if nargin<3
 tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
% Determining the lower bounds. 
S=1:N;
%k=1:n;
%vi=v(bitset(0,k));
%lb=vi';
lb=-inf(n,1);
% Determining the pre-nucleolus
try
   x=msk_prenucl(v);
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
ub=inf(n,1);

prob.a=A2;
prob.blx=lb;
prob.bux=ub;

% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
%param.MSK_IPAR_INTPNT_BASIS   = sc.MSK_OFF;
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex. MSK 7
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8


for ii=1:lpt
    B1=[-v';v(N)];
    C=cmat(ii,:);
    prob.buc=B1;
    prob.c=C;
    [rcode,res] = mosekopt('minimize echo(0)',prob,param);
    sol=res.sol;
    if strcmp(sol.bas.solsta,'OPTIMAL') ~= 1
       sol.bas.solsta
       error('Optimization terminated not successfully')
    end
    alpv(ii)=sol.bas.pobjval;
    v(pt(ii))=alpv(ii);
end

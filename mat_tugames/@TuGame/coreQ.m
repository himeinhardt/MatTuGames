function [crq, x, lS]=coreQ(clv,tol)
% COREQ checks the existence of the core of game v.
% 
%
% Usage: [crq x]=coreQ(clv,tol)
%
% Define variables:
%  output:
%  crq      -- Returns 1 (true) or 0 (false).
%  x        -- The smallest allocation that satisfies all core constraints.
%  lS       -- Weights (dual solution)
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
%   10/29/2012        0.3             hme
%   05/14/2015        0.7             hme
%                


if nargin<2
 tol=10^7*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
S=1:N;
PlyMat=zeros(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k)==1;end
v1(S)=v(S)-tol;
v1(N)=v(N);

% Options
opts.Display='off';
opts.Simplex='on';
opts.LargeScale='on';
opts.Algorithm='dual-simplex';
opts.TolFun=1e-10;
opts.TolX=1e-10;
opts.TolRLPFun=1e-10;
%% for dual-simplex
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);

[x,fval,exitflag,~,lambda]=linprog(ones(n,1),-PlyMat,-v1,[],[],[],[],[],opts);
x=x';
lS=lambda.ineqlin';
crq=abs(v(N)-fval)<=abs(tol);

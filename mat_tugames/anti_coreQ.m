function [acrq x]=anti_coreQ(v,tol)
% ANTI_COREQ checks the existence of the anti-core of game v.
% 
%
% Usage: [acrq x]=anti_coreQ(v,tol)
% Define variables:
%  output:
%  acrq     -- Returns 1 (true) or 0 (false).
%  x        -- The smallest allocation that satisfies all anti-core constraints.
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
%   07/10/2010        0.1 alpha       hme
%   07/09/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%                



if nargin<2
 tol=10^5*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
S(N)=[];
v1(S)=v(S)+tol;
v1(N)=v(N);

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


%opts=optimset('Simplex','on','LargeScale','on','Algorithm','simplex','MaxIter',128);
[x,fval,exitflag]=linprog(-ones(n,1),PlyMat,v1,[],[],[],[],[],opts);
x=x';
%exitflag
%fval
acrq=abs(v(N)-fval)<=abs(tol);


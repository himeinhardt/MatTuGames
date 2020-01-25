function [crq,x,lS]=coreQ(v,tol)
% COREQ checks the existence of the core of game v.
% 
%
% Usage: [crq x]=coreQ(v,tol)
% Define variables:
%  output:
%  crq      -- Returns 1 (true) or 0 (false).
%  x        -- The smallest allocation that satisfies all core constraints.
%  lS       -- Weights (dual solution)
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
%   05/12/2015        0.7             hme
%                


if nargin<2
 tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N;
PlyMat=zeros(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k)==1;end
v1(S)=v(S)-tol;
v1(N)=v(N);
%A=-PlyMat;
A=sparse(-PlyMat);
B=-v1;

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


[x,fval,exitflag,~,lambda]=linprog(ones(n,1),A,B,[],[],[],[],[],opts);
lS=lambda.ineqlin';
x=x';
crq=abs(v(N)-fval)<=abs(tol);

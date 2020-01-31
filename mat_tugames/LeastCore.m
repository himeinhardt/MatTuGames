function [fmin,x,A1,B1,exitflag]=LeastCore(v,tol)
% LEASTCORE computes the least core of game using 
% the Matlab's Optimization toolbox. Uses now Dual-Simplex (Matlab R2015a).
%
%
% Usage: [epsv, x]=LeastCore(v,tol)
% Define variables:
%  output:
%  fmin         -- Critical value of the least core.
%  x            -- An allocation that satisfies the least-core constraints.
%  A1           -- least-core matrix.
%  B1           -- boundary vector.
%  exitflag     -- Exit flag of the LP. Returns 1 for a successful termination.
%
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
%   12/22/2012        0.3             hme
%   09/22/2014        0.5             hme
%   03/28/2015        0.7             hme
%   04/21/2019        1.0             hme
%                


if nargin<2
 tol=10^8*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N;
A1=zeros(N,n);
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[-v';v(N)];
C=[zeros(1,n),1];
ub=[ones(1,n)*v(N),Inf];
%opts=optimset('Simplex','on','LargeScale','off','MaxIter',256);
opts.Display='off';
opts.Simplex='on';
%opts.ActiveSet='on';
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

[xmin,fmin,exitflag,output,lambda]=linprog(C,A2,B1,[],[],[],ub,[],opts);
x=xmin';
x(end)=[];

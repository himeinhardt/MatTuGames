function [crq x]=p_coreQ(v,tol)
% P_COREQ checks the existence of the core of game v using Matlab's PCT..
% 
%
% Usage: [crq x]=p_coreQ(v,tol)
% Define variables:
%  output:
%  crq      -- Returns 1 (true) or 0 (false).
%  x        -- The smallest allocation that satisfies all core constraints.
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
%   11/11/2012        0.3             hme
%                



if nargin<2
 tol=10^8*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
S(N)=[];
v1(S)=v(S)-tol;
v1(N)=v(N);

spmd
 [x,fval,exitflag]=linprog(ones(n,1),-PlyMat,-v1);
end
x=x{1}';
crq=abs(v(N)-fval{1})<=abs(tol);


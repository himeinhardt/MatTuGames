function [acrq a_x]=CddAntiCoreQ(v,tol)
% CDDANTICOREQ checks the existence of the anti-core of game v.
% 
%
% Usage: [acrq a_x]=CddAntiCoreQ(v,tol)
% Define variables:
%  output:
%  crq      -- Returns 1 (true) or 0 (false).
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
%   11/24/2012        0.3             hme
%                



if nargin<2
 tol=10^8*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=PlyMat(1:N,:);
A1(N+1,:)=-PlyMat(end,:);
B1=[v';-v(N)];
objective=PlyMat(end,:);
IN=struct('obj',objective,'A',A1,'B',B1);
OUT = cddmex('solve_lp_DS',IN);
if OUT.how~=1
  acrq = false;
else
  a_x=OUT.xopt;
  acrq=abs(v(N)-OUT.objlp)<=abs(tol);
end

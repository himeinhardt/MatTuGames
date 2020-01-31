function [crq, x]=CddCoreQ(clv,tol)
% CDDCOREQ checks the existence of the core of game v using cddmex.
% 
%
% Usage: [crqi, x, lS]=clv.CddCoreQ(tol)
% Define variables:
%  output:
%  crq      -- Returns 1 (true) or 0 (false).
%  x        -- The smallest allocation that satisfies all core constraints.
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
%   11/22/2012        0.3             hme
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
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat(1:N,:);
A1(N+1,:)=PlyMat(end,:);
B1=[-v';v(N)];
objective=PlyMat(end,:);
IN=struct('obj',objective,'A',A1,'B',B1);
OUT = cddmex('solve_lp_DS',IN);
% lS=-OUT.lambda';
if OUT.how~=1
  crq = false;
else
  x=OUT.xopt';
  crq=abs(v(N)-OUT.objlp)<=abs(tol);
end


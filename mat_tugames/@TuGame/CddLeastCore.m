function [epsv, x1, A1, B1]=CddLeastCore(clv,tol)
% CDDLEASTCORE computes the least core of game clv using cddmex.
% 
%
% Usage: [epsv, x]=CddLeastCore(v,tol)
% Define variables:
%  output:
%  epsv     -- Critical value of the least core.
%  x1       -- An allocation that satisfies the least-core constraints.
%  A1       -- least-core matrix.
%  B1       -- boundary vector.
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
%   12/22/2012        0.3             hme
%                



if nargin<2
 tol=10^8*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat(1:N,:);
A1(N+1,:)=PlyMat(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
B1=[-v';v(N)];
objective=[zeros(1,n),1];
IN=struct('obj',objective,'A',A1,'B',B1);
OUT = cddmex('solve_lp_DS',IN);
x=OUT.xopt;
x1=x';
x1(end)=[];
epsv=OUT.objlp;

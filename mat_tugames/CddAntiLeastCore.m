function [epsv, x1, A1, B1]=CddAntiLeastCore(v,tol)
% CDDANTILEASTCORE computes the least core of game v using cddmex.
%
%
% Usage: [epsv, x]=CddAntiLeastCore(v,tol)
% Define variables:
%  output:
%  epsv     -- Critical value of the anti least core.
%  x        -- An allocation that satisfies the anti least-core constraints.
%  A1       -- anti least-core matrix.
%  B1       -- boundary vector.
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
%   08/09/2016        0.9             hme
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
A1(:,end+1)=1;
A1(N:N+1,end)=0;
B1=[v';-v(N)];
objective=[zeros(1,n),-1];
IN=struct('obj',objective,'A',A1,'B',B1);
OUT = cddmex('solve_lp_DS',IN);
x=OUT.xopt;
x1=x';
x1(end)=[];
epsv=OUT.objlp;

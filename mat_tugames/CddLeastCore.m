function [epsv, x1, A1, B1]=CddLeastCore(v,tol)
% CDDLEASTCORE computes the least core of game v using cddmex.
%
%
% Usage: [epsv, x]=CddLeastCore(v,tol)
% Define variables:
%  output:
%  epsv     -- Critical value of the least core.
%  x        -- An allocation that satisfies the least-core constraints.
%  A1       -- least-core matrix.
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
%   12/12/2012        0.3             hme
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
B1=[-v';v(N)];
objective=[zeros(1,n),1];
IN=struct('obj',objective,'A',A1,'B',B1);
OUT = cddmex('solve_lp_DS',IN);
x=OUT.xopt;
x1=x';
x1(end)=[];
epsv=OUT.objlp;

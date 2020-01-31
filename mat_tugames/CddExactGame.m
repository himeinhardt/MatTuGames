function [v_e, xm]=CddExactGame(v,tol)
% CDDEXACTGAME computes the exact game from v using cddmex.
%
%  Usage: [v_e xm status]=CddExactGame(v,tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
%    
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/25/2016        0.8             hme
%                



if nargin<2
 tol=10^8*eps; % No effect here!
end

if CddCoreQ(v,tol)==0
  error('Core is empty!');
 else
end

N=length(v);
[~, n]=log2(N);
S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat;
A1(N+1,:)=PlyMat(end,:);
B1=[-v';v(N)];
v_e=zeros(1,N);
xm=zeros(N,n);
for S=1:N
    objective=PlyMat(S,:);;
    IN=struct('obj',objective,'A',A1,'B',B1);
    OUT = cddmex('solve_lp_DS',IN);
    x=OUT.xopt;
    xm(S,:)=x';
    v_e(S)=OUT.objlp;
end

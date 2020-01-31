function [v_e, xm]=CddExactGame(clv,tol)
% CDDEXACTGAME computes the exact game from v using cddmex.
%
%  Usage: [v_e xm status]=clv.CddExactGame(tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
%    
%  input:
%  clv      -- TuGame class object.
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

if clv.CddCoreQ(tol)==0
  error('Core is empty!');
 else
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

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

function [v_e xm]=p_glpk_exact_game(v,tol);
% P_EXACT_GAME computes the exact game from game v using glpkmex and 
% Matlab's PCT.
%
% http://www.gnu.org/software/glpk/glpk.html 
%
%  Usage: [v_e xm status]=p_glpk_exact_game(v,tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
%  exfl     -- Returns 5 if the optimization terminated successfully. 
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to -10^8*eps.


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
 tol=10^6*eps;
end

N=length(v);
[l1, n]=log2(N);

if CddCoreQ(v,tol)==0
  error('Core is empty!');
 else
end


s0='U';
ctype=repmat(s0,1,N+1);
ctype(N:N+1)='S';
%ctype=[];
vartype=[];
s=1; % minimization problem 
param.msglev = 1;
param.lpsolver = 1; % Simplex method.



S=1:N;
parfor k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat(1:N,:);
A1(N+1,:)=PlyMat(end,:);
A2=sparse(A1);
B1=[-v';v(N)];
ub=inf(1,n);
v_e=zeros(1,N);
xm=zeros(N,n);
parfor S=1:N
    C=PlyMat(S,:);
    [xmin,fmin,status]=glpk(C,A2,B1,[],ub,ctype,vartype,s,param);
    xm(S,:)=xmin';
    v_e(S)=fmin;
    if status~=5 % old 180
       warning('Solution not optimal!')
    end
end    



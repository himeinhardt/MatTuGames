function [v_e xm]=msk_exact_game(v,tol);
% MSK_EXACT_GAME computes the exact game from game v using mosekmex.
% 
% MSK-SOLVER: http://www.mosek.com/
% 
%  Usage: [v_e xm]=msk_exact_game(v,tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
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
[~,n]=log2(N);

if CddCoreQ(v,tol)==0
  error('Core is empty!');
 else
end

k=1:n;
vi=v(bitset(0,k));

S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat;
A1(N+1,:)=PlyMat(end,:);
A2=sparse(A1);
prob.a=A2;
B1=[-v';v(N)];
ub=inf(n,1);
lb=vi';
v_e=zeros(1,N);
xm=zeros(N,n);
prob.buc=sparse(B1);
prob.blc=-inf(N+1,1);
prob.blx=lb;
prob.bux=ub;
% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
param.MSK_IPAR_INTPNT_BASIS   = 'MSK_ON';
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex. MSK 7
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
%param.MSK_DPAR_BASIS_TOL_X = 1.0e-9;
%param.MSK_DPAR_BASIS_TOL_S = 1.0e-9;
%param=[];

for S=1:N
    C=PlyMat(S,:);
    prob.c=C';
    [rcode,res] = mosekopt('minimize echo(0)',prob,param);
    sol=res.sol;
    x=sol.bas.xx';
    xm(S,:)=x';
%    v_e(S)=C*x';
    v_e(S)=res.sol.bas.pobjval;
end    



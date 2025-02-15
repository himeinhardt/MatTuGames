function [v,x]=msk_linear_production(A,mB,p)
% MSK_LINEAR_PRODUCTION computes from a production problem (A,mB,p) a linear production game using mosekmex.
%
% SOURCE: G. Owen, "On the core of linear production games", Mathematics of
%            Operations Research vol. 9, (1975) 358{370.
%
% Usage: [v,x]=msk_linear_production(A,mB,p)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- A matrix of optimal solution vectors for all
%              coalition S.
% input: 
%  A        -- Production coefficient matrix of size (qxr) w.r.t. inequalities.
%  mB       -- Resource matrix of size (nxq) w.r.t. inequalities.
%  p        -- Price vector of produced goods of size (1xr).
%
%
% Example:
% Define a production coefficient matrix A by
% A =
%    0.2500   0.1000
%    8.0000   5.0000
%    4.0000   6.0000
%
% and a resource matrix mB by
% mB = 
%     9   260   200
%     5   120   200
%    14   590   400    
%
% Finally, consider the price vector of the r-goods through
% p = 
%     200   150    
%
% Now, invoke
%  v=msk_linear_production(A,mB,p)
% to get
% v =
%    1.0e+04 *
%
%    0.7000    0.3600    1.1000    1.4000    2.2091    1.9500    2.6500
%
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/17/2021        1.9             hme
%


[q,r]=size(A);
[n,~]=size(mB);
N=2^n-1;
v=zeros(1,N);
x=zeros(N,r);
l=1:q;

% objective function
prob.c=p';  % max problem
% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
%param.MSK_IPAR_INTPNT_BASIS   = sc.MSK_OFF;
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex. MSK 7
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
%param.MSK_DPAR_BASIS_TOL_X = 1.0e-9;
%param.MSK_DPAR_BASIS_TOL_S = 1.0e-9;

%% boundary on variables.
prob.blx=zeros(r,1);
prob.bux=inf(r,1);

% Setting constraints;
adm=zeros(q,N);
for kk=1:q
    qrs=mB(:,kk)';
    adm(kk,:)=additive_game(qrs);
end



% Sub-Problems
it=0:-1:1-n;
k=1:n;
for S=1:N
   % Matrix and boundary vector for coalition S.
   B=adm(:,S)';
   prob.a=A;
   prob.buc=B;
   prob.blc=-inf(q,1);
   [rcode,res] = mosekopt('maximize echo(0)',prob,param);
   sol=res.sol;
   xx=sol.bas.xx;
   v(S)=p*xx;
   x(S,:)=xx';
end

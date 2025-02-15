function [v,x]=msk_linprog_game(A,H,mB,mD,p)
% MSK_LINPROG_GAME computes from a production matrix (A,H,mB,mD,p) a linear programming game using mosekmex.
%
% SOURCE: G. Owen, "On the core of linear production games", Mathematics of
%            Operations Research vol. 9, (1975) 358{370.
%
% Usage: [v,x]=msk_linprog_game(A,H,mB,md,p)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- A matrix of optimal solution vectors for all
%              coalition S.
% input: 
%  A        -- Production coefficient matrix of size (qxr) w.r.t. inequalites.
%  H        -- Production coefficient matrix of size (q1xr) w.r.t. equalities.    
%  mB       -- Resource matrix of size (nxq) w.r.t. inequalites.
%  mD       -- Resource matrix of size (nxq1) w.r.t. equalities.    
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
%  v=msk_linprog_game(A,[],mB,[],p)
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

if isempty(A)==0 && isempty(H)==0
   [q,r]=size(A)
   [q1,~]=size(H);
   [n,~]=size(mB);
elseif isempty(A)==0 && isempty(H)==1
   [q,r]=size(A);
   [n,~]=size(mB);    
   q1=0;
elseif isempty(A)==1 && isempty(H)==0
   [q,r]=size(H);    
   [n,~]=size(mD);
   q1=0;
else
   error("At least one Matrix must be defined!");
   return;
end    
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

if isempty(H)==0 && isempty(A) ==0
  ctypes(l+1:q+q1)=s1
  ddm=zeros(q1,N);
  for kk=1:q1
      qrs2=-mD(:,kk)';    
      ddm(kk,:)=additive_game(qrs2);
  end
elseif isempty(H)==0 && isempty(A) ==1
  ctype(l)=s1;  
  ddm=zeros(q,N);
  for kk=1:q
      qrs2=-mD(:,kk)';    
      ddm(kk,:)=additive_game(qrs2);
  end    
else
   ddm=[]; 
end

% Setting constraints;
if isempty(A)==0
  adm=zeros(q,N);
   for kk=1:q
       qrs=mB(:,kk)';
       adm(kk,:)=additive_game(qrs);
   end
else
   adm=[]; 
end    

% Sub-Problems
it=0:-1:1-n;
k=1:n;
if isempty(H)==1 && isempty(A)==0
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
elseif isempty(H)==0 && isempty(A)==1
  for S=1:N
       % Matrix and boundary vector for coalition S.
       D=ddm(:,S)';
       prob.a=H;
       prob.buc=D;
       prob.blc=D;
       [rcode,res] = mosekopt('maximize echo(0)',prob,param);
       sol=res.sol;
       xx=sol.bas.xx;
       v(S)=p*xx;
       x(S,:)=xx';
  end    
elseif isempty(H)==1 && isempty(A)==1   
  for S=1:N
       % Matrix and boundary vector for coalition S.
       A1=[A;H];
       B=adm(:,S)';
       D=ddm(:,S)';
       prob.a=A1;
       prob.buc=[B;D];
       prob.blc=[-inf(q1,1);D];
       [rcode,res] = mosekopt('maximize echo(0)',prob,param);
       sol=res.sol;
       xx=sol.bas.xx;
       v(S)=p*xx;
       x(S,:)=xx';
  end
end  

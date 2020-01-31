function [v xm]= msk_AssignmentGame(slv,pfm)
% MSK_ASSIGNMENTGAME computes from an assignment problem (sl_vec,prof_mat) 
% the corresponding symmetric assignment game using mosekmex. If the problem is 
% not symmetric, it will be transformed into a symmetric one.
% The assignment game will be derived from an integer problem. 
%
% MSK-SOLVER: http://www.mosek.com/
%
%
% Usage: [v xm]=msk_AssignmentGame(sl_vec,prof_mat)
%
% Define variables:
%  output:
%  v        -- An assignment game v of length 2^n-1.
%  xm       -- Assignment matrix. 
%              (optimal solution vectors for each integer problem).
%            
%  input:
%  sl_vec   -- A vector of sellers. From this vector, the vector of 
%              buyers will be constructed.
%              Example: sl_vec=[1 2 3 4 5];
%  prof_mat -- A square profit matrix.
%              use the function profit_matrix() to compute one.
%              or invoke
%              prfm=magic(5);
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/11/2014        0.5             hme
%
nv=length(slv);
n=2*nv;
sq=nv^2;
N=2^n-1;
S=1:N;
pl=1:n;
%
% Create integer problem
trp=pfm';
f=trp(:)';
idm=eye(nv);
A1=zeros(n,sq);
n1=1;
en=nv;
for k=1:nv
   A1(k,n1:en)=1;
   A1(nv+1:end,n1:en)=idm;
   n1=en+1;
   en=en+nv;
end
v=zeros(1,N);
A=sparse(A1);
clear A1;

prob.a=A;
prob.c=-f;
%prob.blc=-inf(N+1,1);
prob.blx=zeros(sq,1);
prob.bux=ones(sq,1);
%prob.bux=ub;
[rcode,res] = mosekopt('param echo(0)');
param=res.param;

if nargout > 1
   xm=zeros(N,sq);
end

% Integer problem for each coalition S.
for jj=1:N
   b=bitget(jj,pl);
   % [x fval exitflag] = bintprog1(-f,A,b);
   prob.buc=b;
   [rcode,res] = mosekopt('minimize echo(0)',prob,param);
   sol=res.sol;
   x=sol.bas.xx;
   v(jj)=f*x;
   if nargout > 1
       xm(jj,:)=x';
   end
end

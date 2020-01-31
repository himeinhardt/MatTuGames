function [v xm]= p_msk_AssignmentGame(slv,pfm)
% P_MSK_ASSIGNMENTGAME computes from an assignment problem (sl_vec,prof_mat) 
% the corresponding symmetric assignment game using mosekmex and Matlab's PCT.
% If the problem is  not symmetric, it will be transformed into a symmetric one.
% The assignment game will be derived from an integer problem. 
%
% MSK-SOLVER: http://www.mosek.com/
%
%
% Usage: [v xm]=p_msk_AssignmentGame(sl_vec,prof_mat)
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
%   05/16/2019        1.1             hme
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

for kk=1:N
  Prob(kk).a=A;
  Prob(kk).c=-f;
  Prob(kk).blx=zeros(sq,1);
  Prob(kk).bux=ones(sq,1);
  Prob(kk).buc=bitget(kk,pl);;
end

ng=nargout;

if ng > 1
   xm=zeros(N,sq);
end


% Integer problem for each coalition S.
parfor jj=1:N
   prob=Prob(jj);
   [rcode,res] = mosekopt('param echo(0)');
   param=res.param;
   param.MSK_IPAR_NUM_THREADS = 1;
   [rcode,res] = mosekopt('minimize echo(0)',prob,param);
   sol=res.sol;
   x=sol.bas.xx;
   v(jj)=f*x;
   if ng > 1
       xm(jj,:)=x';
   end
end

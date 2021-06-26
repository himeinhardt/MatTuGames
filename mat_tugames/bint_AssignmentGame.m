function [v xm exfl]= bint_AssignmentGame(slv,pfm)
% BINT_ASSIGNMENTGAME computes from an assignment problem (sl_vec,prof_mat) 
% the corresponding symmetric assignment game. If the problem is 
% not symmetric, it will be transformed into a symmetric one.
% The assignment game will be derived from an integer problem.
%
% Usage: [v xm exfl]=bint_AssignmentGame(sl_vec,prof_mat)
%
% Define variables:
%  output:
%  v        -- An assignment game v of length 2^n-1.
%  xm       -- Assignment matrix. 
%              (optimal solution vectors for each integer problem).
%  exfl     -- Vector of exitflags.
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
%   10/21/2014        0.5             hme
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
intcon=1:sq;
A=sparse(A1);
clear A1;
lb=zeros(sq,1);
ub=ones(sq,1);
Aeq=[];
beq=[];
b=[];
opts=optimoptions('intlinprog','Display','off');
%opts.BranchingRule = 'mostfractional';
opts.CutGeneration = 'intermediate';
opts.CutGenMaxIter = 20; % default is 10
% opts.CutGeneration = 'none';
opts.TolInteger = 1e-06;
opts.TolGapAbs = 0;

if nargout > 1
   xm=zeros(N,sq);
   exfl=false(1,N);
end

% Define binary integer problem
    problem = struct('f',-f,'intcon',intcon,...
    'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,...
    'lb',lb,'ub',ub,'options',opts,...
    'solver','intlinprog');

% Integer problem for each coalition S.
for jj=1:N
   b=bitget(jj,pl);
   problem.bineq=b;
   [x fval exitflag output] = intlinprog(problem);
   v(jj)=-fval;
   if nargout > 1
       xm(jj,:)=x';
       exfl(jj)=exitflag;
   end
end

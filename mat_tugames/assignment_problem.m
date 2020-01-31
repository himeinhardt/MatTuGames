function problem = assignment_problem(slv,pfm,ss)
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
%  ss       -- A coalition S represented by its unique integer.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/10/2014        0.5             hme
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
f=trp(:)'
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
intcon=pl;
A=sparse(A1);
%clear A1;
lb=zeros(sq,1);
ub=ones(sq,1);
Aeq=[];
beq=[];
opts=optimoptions('intlinprog','Display','final');
%opts.BranchingRule = 'mostfractional';
%opts.CutGeneration = 'intermediate';
%opts.CutGenMaxIter = 20; % default is 10
% opts.CutGeneration = 'none';
%opts.TolInteger = 1e-06;
%opts.TolGapAbs = 0;


% Integer problem for each coalition S.
   b=bitget(ss,pl);
   %[x fval exitflag] = bintprog1(-f,A,b);
   % Problem unbounded with intlinprog MATLAB R2014a
   %try
   problem = struct('f',-f,'intcon',pl,...
    'Aineq',A1,'bineq',b,'Aeq',Aeq,'beq',beq,...
    'lb',lb,'ub',ub,'options',opts,...
    'solver','intlinprog');
   %[x fval exitflag output] = intlinprog(problem);
   %fval
   %output
   % [x fval exitflag] = intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub,opts);
   %[x fval exitflag] = intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub);
   %catch 
   %  [x fval exitflag] = bintprog(-f,A,b);
   %end

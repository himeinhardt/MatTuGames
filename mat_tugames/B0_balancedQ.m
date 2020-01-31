function [bcQ, cmat, rk, cf]=B0_balancedQ(cS,n,b0,tol)
% B0_BALANCEDQ verifies whether the collection of coalitions is weakly balanced. 
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [bcQ, cmat, rk]=B0_balancedQ(cS,n,b0,tol)
%
% 1. Example:
% Choose a collection of coalitions:
%
% cB={[2],[1 3],[1 4],[2 3 4]}
% uB=clToMatlab(cB)                   
% 
%  uB =
%
%     2     5     9    14
%
% [bcQ, ~, ~,cf]=B0_balancedQ(uB,4)
%
% bcQ =
%
%     1
%
% cf =
%
%    0.5000
%    0.5000
%    0.5000
%    0.5000
%
%
% 2. Example:
% A collection of sets given by their unique integer representation:
%
% cS=[1   254   253     2    16   239   252     3   127   191   223   247   251     4];
% n=8;
% bSQ=B0_balancedQ(cS,n) 
% 
% bSQ =
%
%     1
% 
% Define variables:
%  output:
%  bcQ      -- Returns 1 (true) or 0 (false).
%  cmat     -- Incidence matrix of players. 
%  rk       -- Rank of matrix cmat.
%  cf       -- Balanced weights.
%
%
%  input:
%  cS        -- Collection of coalitions.
%  n         -- Number of players involved.
%  b0        -- Collection of coalitions of size 1, can be empty. However,
%               it is recommended to supply this set. In this case use
%               as input a cell like b0={[1],[4]}, which is b0=[1,8]
%               in the unique integer representation. If one supplies
%               [1 4] instead one will get wrong results.
%  tol       -- Tolerance value. Its default value is set to 10^4*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/03/2017        0.9             hme
%   02/24/2019        1.0             hme
%                

    
if nargin < 3
   b0=[]; 
   tol=10^4*eps;
elseif nargin < 4
   tol=10^4*eps;
end

bcQ=false;

N=2^n-1;

zv=zeros(n,1);
if isempty(b0)
   k=1:n;
   ic=2.^(k-1);
else
   if iscell(b0)
      ic=clToMatlab(b0);
   else
      ic=b0;
   end
end
pws=PowerSet(ic);
lic=length(pws);
ls=lic+1;
sS=cell(1,ls);
for k=1:lic
   sS{k}=unique([cS,pws{k}]);
end
sS{ls}=cS;


warning('off','all');
for k=1:ls
%%sS{k}
   if nargout < 4
      [cmat,xS,ef]=CheckB0Bal(n,sS{k},tol);
   elseif nargout == 4
   % Trying to find positive weights.
      [cmat,xS,ef,cf]=CheckB0Bal(n,sS{k},tol); 
   else
      [cmat,xS,ef]=CheckB0Bal(n,sS{k},tol);
   end
  if ef~=1
  else
      bcQ=all(abs(zv-xS)<tol);
  end
  rk=rank(cmat);
  if bcQ==1 
     warning('on','all');
     return;
  end
end
warning('on','all');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cmat,sol,ef,cf]=CheckB0Bal(n,iS,tol)
% CHECKBAL checks balancedness of the collection iS.
%
%
% Define variables:
%  output:
%  cmat      -- Incidence matrix of players.
%  sol       -- Solution vector (zeros). 
%  ef        -- Exitflag of the linear problem.
%  cf        -- Balanced weights.
%
%  input:
%  n         -- Number of players involved.
%  iS        -- Collection of coalitions.
%  tol       -- Tolerance value. Its default value is set to 10^4*eps.
%

int=0:-1:1-n;
ov=ones(n,1);
N=2^n-1;

liS=length(iS);
cmat=(rem(floor(iS(:)*pow2(int)),2)==1)';
cmat=double(cmat);
[c1,c2]=size(cmat);
A=-cmat';
ovn=ones(c2,1);
b=zeros(c2,1);
zf=A'*ovn;
f=zf';
Aeq=ov';
beq=0;
mtv=verLessThan('matlab','9.1.0');
    try
      if mtv==1
         options = cplexoptimset('MaxIter',128,'Dual-Simplex','on','Display','off');
      else
         options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
      end
      options.barrier.convergetol=1e-12;
      options.simplex.tolerances.feasibility=1e-9;
      options.simplex.tolerances.optimality=1e-9;
      options.emphasis.numerical=1;
      options.barrier.display=0;
      options.feasopt.tolerance=1e-12;
      options.Param.lpmethod=2;
      [sol,fval,ef,~,lambda] = cplexlp(f,A,b,Aeq,beq,[],[],[],options);
    catch
      opts.Display='off';
      opts.Simplex='on';
      opts.LargeScale='on';
      opts.Algorithm='dual-simplex';
      opts.TolFun=1e-10;
      opts.TolX=1e-10;
      opts.TolRLPFun=1e-10;
      %% for dual-simplex
      opts.MaxTime=9000;
      opts.Preprocess='none';
      opts.TolCon=1e-6;
      opts.MaxIter=10*(N+n);
      [sol,fval,ef,~,lambda] = linprog(f,A,b,Aeq,beq,[],[],[],opts);
    end
%
% Trying to find positive weights.
%
if ef==1
  y1=lambda.ineqlin;
  y2=lambda.eqlin;
  cf=(y1+1);
%% Cross check result must give the zero vector.
  bq=-cmat*cf+y2*ov;
else
  cf=[];
end

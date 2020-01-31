function [bcQ, cmat, S]=p_AntiBalancedSetQ(exc,n,S,tol)
% P_ANTIBALANCEDSETQ verifies whether the set of induced
% coalitions is a anti balanced collection. Checking Kohlberg's criterion.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a). It is recommended not to 
% use this function, use p__B0_balancedCollectionQ instead.
% 
%
% Usage: bcQ=p_AntiBalancedSetQ(exc,n,S,tol)
% Define variables:
%  output:
%  bcQ      -- Returns 1 (true) or 0 (false).
%  cmat     -- Incidence matrix of players.
%  S        -- Collection of coalitions (unique integer representation)
%
%  input:
%  exc      -- The excess of game v w.r.t. x of length 2^n-2. 
%  n        -- Size of underlying player set.
%  S        -- A collection of sets. 
%  tol      -- Tolerance value. Its default value is set to 10^4*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/06/2017        0.9             hme
%                

    
if nargin < 3 
   tol=10^4*eps;
elseif nargin < 4
   tol=10^4*eps;
end

bcQ=false;

zv=zeros(n,1);
[sx,idx]=sort(exc,'descend');
mx=min(sx);
slc=sx<=mx+tol;
% S=idx(slc);
thex=mx;
rsx=sx;
rsx(slc)=[];
idx(slc)=[];

warning('off','all');
while 1
   [cmat,xS,ef]=p_CheckAntiBal(n,S);
   if isempty(S)
      bcQ=false;
      warning('on','all');
      break;
   elseif ef~=1
      bcQ=false;
      warning('on','all');
      break;
   end
    bcQ=all(abs(zv-xS)<tol);
    rk=rank(cmat);
   if bcQ==0
      warning('on','all');
      break;
   elseif rk==n && bcQ==1
      warning('on','all');
      break;
   end
   slc=rsx<=thex;
   sS=idx(slc);
   idx(slc)=[];
   rsx(slc)=[];
   S=[sS,S];
   mx=min(rsx);
   thex=mx+tol;
end


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cmat,sol,ef]=p_CheckAntiBal(n,iS)
% CHECKANTIBAL checks balancedness of the collection iS.
%
%
% Define variables:
%  output:
%  cmat      -- Incidence matrix of players.
%  sol       -- Solution vector (zeros).
%  ef        -- Exitflag of the linear problem.
%
%  input:
%  n         -- Number of players involved.
%  iS        -- Collection of coalitions.
%

int=0:-1:1-n;
ov=ones(n,1);
N=2^n-1;

liS=length(iS);
cmat=(rem(floor(iS(:)*pow2(int)),2)==1)';
cmat=double(cmat);
[c1,c2]=size(cmat);
A=-cmat';
%rk=rank(A)
ovn=ones(c2,1);
b=zeros(c2,1);
zf=A'*ovn;
f=zf';
Aeq=ov';
beq=0;
    try 
      options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
      options.barrier.convergetol=1e-12;
      options.simplex.tolerances.feasibility=1e-9;
      options.simplex.tolerances.optimality=1e-9;
%      options.threads=16;
      options.emphasis.numerical=1;
      options.barrier.display=0;
      options.feasopt.tolerance=1e-12;
      options.Param.lpmethod=2;
      [sol,fval,exitflag,~,lambda] = cplexlp(f,A,b,Aeq,beq,[],[],[],options);
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
      [sol,fval,exitflag,~,lambda] = linprog(f,A,b,Aeq,beq,[],[],[],opts);
    end
ef=exitflag; 



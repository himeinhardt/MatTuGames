function [bcQ, cmat, iS]=weightedCsBalancedCollectionQ(clv,cs,x,pS,tol)
% WEIGHTEDCSBALANCEDCOLLECTIONQ verifies whether the set of induced
% coalitions is a weighted balanced collection w.r.t. the coaltion structure cs. 
% Checking Kohlberg's criterion.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
% 
%
% Usage: bcQ=clv.weightedCsBalancedCollectionQ(cs,x,pS,tol)
% Define variables:
%  output:
%  bcQ      -- Returns 1 (true) or 0 (false).
%  cmat     -- Incidence matrix of players.
%  iS       -- Collection of coalitions (unique integer representation)
%
%  input:
%  clv      -- TuGame class object.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  x        -- A (pre-)imputation of length n.
%  pS       -- A vector of weights of length 2^n-1 (optional).
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/24/2017        0.9             hme
%   05/26/2024        1.9.2           hme
%                

if nargin < 4 
   tol=10^8*eps;
   pS='';
elseif nargin < 5
   tol=10^8*eps;
end

bcQ=false;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
if iscell(cs)
   cs=clToMatlab(cs);
else
  cs=double(cs);
end
exc=clv.excess(x);
ecs=exc(cs);
effQ=all(abs(ecs)<tol);

if effQ==0
   cmat=[];
   S=[];
   return;
end


if isempty(pS)
  int=1-n:1:0;
  S=1:N;
  mat=(rem(floor(S(:)*pow2(int)),2)==1)';
  clS=ones(1,n)*mat;
  pS=1./clS; % weights of coalitions (per capita)
end

zv=zeros(n,1);
exc=clv.weightedExcess(x,pS);
[sx,idx]=sort(exc,'descend');
sli=ismember(idx,cs)==0;
sx=sx(sli);
idx=idx(sli);
mx=max(sx);
slc=sx>=mx-tol;
S=idx(slc);
thex=mx;
rsx=sx;
rsx(slc)=[];
idx(slc)=[];

warning('off','all');
while 1
   [cmat,xS,ef]=cs_WhgsCheckBal(n,S,cs);
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
   slc=rsx>=thex;
   sS=idx(slc);
   idx(slc)=[];
   rsx(slc)=[];
   S=[sS,S];
   mx=max(rsx);
   thex=mx-tol;
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cmat,sol,ef]=cs_WhgsCheckBal(n,iS,cs)
% WHGSCHECKBAL computes the set of induced coalitions.
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
iS=[iS,cs];
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
mtv=verLessThan('matlab','9.1.0');
mtv2=verLessThan('matlab','9.1.12');
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
      [sol,fval,exitflag,~,lambda] = cplexlp(f,A,b,Aeq,beq,[],[],[],options);
    catch
      opts.Display='off';
      opts.Simplex='on';
      opts.LargeScale='on';
      mth1=verLessThan('matlab','24.1.0');
      if mth1==0,
          opts.Algorithm='dual-simplex-highs';
      else
          opts.Algorithm='dual-simplex';
      end      
      opts.TolFun=1e-10;
      opts.TolX=1e-10;
      opts.TolRLPFun=1e-10;
      %% for dual-simplex
      opts.MaxTime=9000;
      opts.Preprocess='none';
      opts.TolCon=1e-6;
      opts.MaxIter=10*(N+n);
      if mtv2==0
         [sol,fval,exitflag,~,lambda] = linprog(f,A,b,Aeq,beq,[],[],opts);
      else %% old api (before R2022a) with initial value.
         [sol,fval,exitflag,~,lambda] = linprog(f,A,b,Aeq,beq,[],[],[],opts);
      end      
    end
ef=exitflag;

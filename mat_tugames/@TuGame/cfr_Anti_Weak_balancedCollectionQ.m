function [bcQ, cmat, S]=cfr_Anti_Weak_balancedCollectionQ(clv,F,x,tol)
% CFR_ANTI_WEAK_BALANCEDCOLLECTIONQ verifies whether the set of induced
% coalitions is a weak_balanced collection with coalition formation restrictions. 
% Checking Kohlberg's criterion.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
% 
%
% Usage: bcQ=clv.cfr_Anti_Weak_balancedCollectionQ(x,F,tol)
% Define variables:
%  output:
%  bcQ      -- Returns 1 (true) or 0 (false).
%  cmat     -- Incidence matrix of players.
%  S        -- Collection of coalitions (unique integer representation)
%
%  input:
%  clv      -- TuGame class object.
%  F        -- For instance, a characterization set for the nucleolus.
%              F must contain the grand coalition N.
%  x        -- A (pre-)imputation of length n.
%  tol      -- Tolerance value. Its default value is set to 10^4*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/22/2017        0.9             hme
%   05/25/2024        1.9.2           hme
%                


if nargin < 2
    error('A collection of sets F is required!');    
if nargin < 3 
   warning('Using default payoff!'):
   x=clv.cfr_Anti_Nucl(F);
   tol=10^4*eps;
elseif nargin < 4
   tol=10^4*eps;
end

bcQ=false;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
ic=clv.tuvi;

effQ=abs((sum(x)-v(N))<tol);

if effQ==0
   cmat=[];
   return;
end    

zv=zeros(n,1);
exc=clv.excess(x);
nic=N-ic;
%if any(exc(nic)>0)
%   bcQ=false;
%   cmat=[];
%   S=[];
%   return;
%else
  iex=exc(nic)==0;
%end

%% Addining the singleton coalitions to F, 
%% since v is definded over F ans si.
F=unique([F,ic]);

b0=ic(iex);
%exc(N)=[];
cS=1:N;
CS=cS(ismember(cS,F)==0);
%F(end)=[];
exc(CS)=[];
[sx,idx]=sort(exc,'ascend');
idx=F(idx);
mx=min(sx);
slc=sx<=mx+tol;
S=idx(slc);
thex=mx;
rsx=sx;
rsx(slc)=[];
idx(slc)=[];

warning('off','all');
while 1
   [cmat,xS,ef]=AntiWeakCheckB0Bal(n,S,b0,N);
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
function [cmat,sol,ef]=AntiWeakCheckB0Bal(n,iS,b0,N)
% CHECKBAL checks B0-balancedness of the collection iS.
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
iS=[iS,N];
liS=length(iS);
cmat=(rem(floor(iS(:)*pow2(int)),2)==1)';
cmat=double(cmat);
[c1,c2]=size(cmat);
A0=-cmat';
%b0
if isempty(b0)==0
   b0m=(rem(floor(b0(:)*pow2(int)),2)==1)';
   b0m=double(b0m)';
   A=[A0;-b0m];
else
   A=A0;
end
[a1,a2]=size(A);
ovn=ones(c2,1);
b=zeros(a1,1);
zf=A0'*ovn;
f=zf';
% Grand coalition
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
      else  %% old api (before R2022a) with initial value.
          [sol,fval,exitflag,~,lambda] = linprog(f,A,b,Aeq,beq,[],[],[],opts);
      end
    end
ef=exitflag; 



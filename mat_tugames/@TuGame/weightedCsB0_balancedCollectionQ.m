function [bcQ, cmat, S]=weightedCsB0_balancedCollectionQ(clv,cs,x,pS,tol)
% WEIGHTEDCSB0_BALANCEDCOLLECTIONQ verifies whether the set of induced
% coalitions is a B0_balanced collection w.r.t. the coaliton structure cs. 
% Checking Kohlberg's criterion.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
% 
%
% Usage: bcQ=clv.weightedCsB0_balancedCollectionQ(cs,x,pS,tol)
% Define variables:
%  output:
%  bcQ      -- Returns 1 (true) or 0 (false).
%  cmat     -- Incidence matrix of players.
%  S        -- Collection of coalitions (unique integer representation)
%
%  input:
%  clv      -- TuGame class object.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  x        -- A (pre-)imputation of length n.
%  pS       -- A vector of weights of length 2^n-1 (optional).
%  tol      -- Tolerance value. Its default value is set to 10^4*eps.


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
   tol=10^4*eps;
   pS='';
elseif nargin < 5
   tol=10^4*eps;
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
ic=clv.tuvi;
if any(exc(ic)>0)
   bcQ=false;
   cmat=[];
   S=[];
   return;
else
  iex=exc(ic)==0;
end
b0=ic(iex);
pws=PowerSet(b0);
lb0=length(pws);
leb=lb0+1;
cS=cell(1,leb);

[sx,idx]=sort(exc,'descend');
sli=ismember(idx,cs)==0;
sx=sx(sli);
idx=idx(sli);
mx=max(sx);
slc=sx>=mx-tol;
S=idx(slc);
for k=1:lb0
   cS{k}=unique([S,pws{k}]);
end
cS{leb}=S;
thex=mx;
rsx=sx;
rsx(slc)=[];
idx(slc)=[];

warning('off','all');
for kk=1:leb
if bcQ==1
   warning('on','all');
   return;
end
   while 1
       [cmat,xS,ef]=cs_CheckB0Bal(n,cS{kk},cs);
       if isempty(cS{kk})
          bcQ=false;
          break;
       elseif ef~=1
          bcQ=false;
          break;
       end
       bcQ=all(abs(zv-xS)<tol);
       rk=rank(cmat);
       if bcQ==0
          break;
       elseif rk==n && bcQ==1
          S=cS{kk};
          break;
       end
       slc=rsx>=thex;
       sS=idx(slc);
       idx(slc)=[];
       rsx(slc)=[];
       S=cS{kk};
       cS{kk}=[sS,S];
       mx=max(rsx);
       thex=mx-tol;
   end
end


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cmat,sol,ef]=cs_CheckB0Bal(n,iS,cs)
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
iS=[iS,cs];
liS=length(iS);
cmat=(rem(floor(iS(:)*pow2(int)),2)==1)';
cmat=double(cmat);
[c1,c2]=size(cmat);
A=-cmat';
ovn=ones(c2,1);
b=zeros(c2,1);
zf=A'*ovn;
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



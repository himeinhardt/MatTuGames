function [bcQ,cmat,S]=balancedCollectionQ(v,x,tol)
% BALANCEDCOLLECTIONQ verifies whether the set of induced
% coalitions is a balanced collection. Checking Kohlberg's criterion.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
% 
%
% Usage: bcQ=balancedCollectionQ(v,x,tol)
% Define variables:
%  output:
%  bcQ      -- Returns 1 (true) or 0 (false).
%  cmat     -- Incidence matrix of players.
%  S        -- Collection of coalitions (unique integer representation)
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- A (pre-)imputation of length n.
%  tol      -- Tolerance value. Its default value is set to 10^4*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/04/2013        0.5             hme
%   01/03/2015        0.6             hme
%   03/28/2015        0.7             hme
%   02/24/2018        0.9             hme
%                

    
if nargin < 2 
   tol=10^5*eps;
elseif nargin < 3
   tol=10^5*eps;
end

bcQ=false;

N=length(v);
[~, n]=log2(N);
effQ=abs(sum(x)-v(N))<tol;

if effQ==0
   cmat=[];
   S=[];
   return;
end    

exc=excess(v,x);
% we need to comment it out for grasping the case N=1.
%exc(N)=[];
[sx,idx]=sort(exc,'descend');
mx=sx(1);
thex=mx-tol;
slc=sx>=thex;
S=idx(slc);
rsx=sx;
rsx(slc)=[];
idx(slc)=[];

warning('off','all');
while 1
   [cmat,xS,ef]=CheckBal(n,S,N);
   if isempty(S)
      bcQ=false;
      warning('on','all');
      break;
   elseif ef~=1
      bcQ=false;
      warning('on','all');
      break;
   end
    bcQ=all(abs(xS)<tol);
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
function [cmat,sol,ef]=CheckBal(n,iS,N)
% CHECKBAL checks balancedness of the collection iS.
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
A=-cmat';
%rk=rank(A)
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
%lambda
%w=lambda.lower
%z=lambda.upper
%y1=lambda.ineqlin
%y2=lambda.eqlin
%y=[y1;y2];
%-cmat
%bq=-cmat*(y1+1)+y2*ov
%z=bq;
%[cmat,ov]*y + w -z - zf
ef=exitflag; 



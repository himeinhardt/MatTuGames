function [bcQ, cmat, S]=B0_balancedCollectionQ(v,x,tol)
% B0_BALANCEDCOLLECTIONQ verifies whether the set of induced
% coalitions is a B0_balanced collection. Checking Kohlberg's criterion.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
% 
%
% Usage: bcQ=B0_balancedCollectionQ(v,x,tol)
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
%   01/24/2017        0.9             hme
%                

    
if nargin < 2 
   tol=10^4*eps;
elseif nargin < 3
   tol=10^4*eps;
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
k=1:n;
ic=2.^(k-1);
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
lb0=numel(pws);
leb=lb0+1;
cS=cell(1,leb);

exc(N)=[];
[sx,idx]=sort(exc,'descend');
mx=sx(1);
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
       [cmat,xS,ef]=CheckB0Bal(n,cS{kk},N);
       if isempty(cS{kk})
          bcQ=false;
          break;
       elseif ef~=1
          bcQ=false;
          break;
       end
       bcQ=all(abs(xS)<tol);
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
function [cmat,sol,ef]=CheckB0Bal(n,iS,N)
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
A=-cmat';
ovn=ones(c2,1);
b=zeros(c2,1);
zf=A'*ovn;
f=zf';
% Grand coalition
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
ef=exitflag; 



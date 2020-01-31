function [bcQ,cmat,S]=modiclusQ(v,x,tol)
% MODICLUSQ verifies whether the set of induced
% coalitions is a bi-balanced collection. Checking adjusted Kohlberg's 
% criterion w.r.t. the modiclus.
%
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
% 
%
% Usage: bcQ=modiclusQ(v,x,tol)
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
%   12/21/2017        0.9             hme
%                

    
if nargin < 2 
   tol=10^4*eps;
elseif nargin < 3
   tol=10^4*eps;
end

bcQ=false;

N=length(v);
[~, n]=log2(N);
dv=dual_game(v);
N1=N+1;
n1=2*n;
N2=2^n1-1;
geQ=all(v<=dv);
effQ=abs(sum(x)-v(N))<tol;

if effQ==0
   cmat=[];
   S=[];
   return;
end    

exc=excess(v,x);
dxc=excess(dv,x);

bxc=zeros(1,N2);
%ii=1;
for k=1:N1;
    for jj=1:N1
        if geQ
          if k>1 && jj >1
             ii=(k-1)+(jj-1)*N1;
             bxc(ii)=exc(k-1)+dxc(jj-1);
           elseif k==1 && jj >1
             ii=N1*(jj-1);
             bxc(ii)=dxc(jj-1);
           elseif k>1 && jj==1
             ii=k-1;
             bxc(ii)=exc(k-1);
           end
        else
          if k>1 && jj >1
             ii=N1*(k-1)+(jj-1);
             bxc(ii)=exc(k-1)+dxc(jj-1);
           elseif k==1 && jj >1
             ii=jj-1;
             bxc(ii)=dxc(jj-1);
           elseif k>1 && jj==1
             ii=N1*(k-1);
             bxc(ii)=exc(k-1);
           end
        end
    end;
end
bxc(N2)=[];
[sx,idx]=sort(bxc,'descend');
clear bxc;
mx=sx(1);
thex=mx-tol;
slc=sx>=thex;
S=idx(slc);
rsx=sx;
rsx(slc)=[];
idx(slc)=[];

warning('off','all');
while 1
   [cmat,xS,ef]=CheckBal(n1,S,N2);
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
iS=[iS,N];
liS=length(iS);
cmat=(rem(floor(iS(:)*pow2(int)),2)==1)';
cmat=double(cmat);
n2=n/2;
ov=ones(n2,1);
cmat=cmat(1:n2,:)+cmat(n2+1:n,:);
cmat(:,end)=ones(1,n2);
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
ef=exitflag; 



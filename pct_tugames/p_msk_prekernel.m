function [x, Lerr, smat, xarr]=p_msk_prekernel(v,x)
% P_MSK_PREKERNEL computes from (v,x) a pre-kernel element using mosekmex
% and Matlab's PCT.
% Source: Meinhardt, 2010.
%
% MSK-Solver: http://www.mosek.com/  
%
% Usage: [x Lerr smat xarr] = p_msk_prekernel(v,x)
%
% Define variables:
%  output:
%  x        -- Kernel element (output)
%  Lerr     -- List of computed function values of hx and h. 
%  smat     -- Matrix of maximum surpluses.
%  xarr     -- History of computed solution at each iteration step.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/29/2013        0.3             hme
%   10/11/2013        0.5             hme
%   02/17/2017        0.9             hme
%   05/03/2019        1.1             hme
%                

if nargin<1
    error('At least the game must be given!');
elseif nargin<2
    N=length(v);
    [~, n]=log2(N);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
    k=1:n;
    if N==1,
      Si=N;
     else
      Si=bitset(N,k,0);
    end
    mv=max(v);
    mnQ=mv>v(N);
    x=(mv-v(Si))/2;
    sx=sum(x);
     if sx>0
        if mnQ==1
           if n < 10
              x(n)=v(N);
           elseif n >= 11 && n <  15
              x=x*mv/sx;
           else
              x=ones(1,n);
           end
        else
           x=x*mv/sx;
        end
      elseif all(abs(x-0)<10^3*eps)==1
       x=(mv-v(bitset(0,k)))/2;
       sx=sum(x);
       mmq=min(x)~=max(x);
       if sx>0 && mmq
          x=select_starting_pt(v);
        else
         x=(v(N)/n)*ones(1,n);
       end
      else
       x=(v(N)/n)*ones(1,n);
     end
    smc=1;
else
    N=length(v);
    [~, n]=log2(N);
    mv=max(v);
    mnQ=mv>v(N);
    if (2^n-1)~=N
       error('Game has not the correct size!');
    end
    smc=1;
end

if islogical(v)
   v=double(v);
end

[x, Lerr, smat, xarr]=computePrk(v,x,smc,0,mnQ);
x=x';
smat=tril(smat,-1)+triu(smat,1);

% Main function to compute a
% pre-kernel element.
%-----------------------------
function [x, Lerr, smat, xarr]=computePrk(v,x,smc,slv,mnQ)
% 
%  output:  -- as above.
%
%  input:
%   v       -- as above.
%   x       -- as above.
%  smc      -- selecting from effc the smallest/largest cardinality (optional).
%              Value must be set to 0 (largest),1 (smallest) or 2 (reset).
%  slv      -- selecting a different linear solver (QR/SVD-decomposition). 
%              Value must be set to 0 or 1.

n=length(x);
N=2^n-1;
x=x';
cnt=0;
if 15<=n 
 CNT=n+2;
else
 CNT=2*(n+1);
end
Lerr=-inf(CNT,2);
xarr=-inf(CNT,n);
m=1+n*(n-1)/2;
upe=true(n);

% Conic Optimization
[rcode,res] = mosekopt('param symbcon echo(0)');
param=res.param;
param.MSK_IPAR_INTPNT_BASIS   = 'MSK_ON';
%param.MSK_IPAR_OPTIMIZER = 1;  % Using interior-point optimizer (MSK 7).
% param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1.0000e-10;
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1.0000e-10;
param.MSK_DPAR_INTPNT_CO_TOL_MU_RED = 1.0000e-10;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0000e-10;
param.MSK_DPAR_INTPNT_CO_TOL_INFEAS = 1.0000e-10;

% Cycling may occur, so that we need an artificial halt
while cnt<CNT  
    cnt=cnt+1;
    [A, smat]=effCoalitions(v,x,smc,cnt);
    upe=tril(upe,-1);
    etr12=A';
    ec12=etr12(upe)';
    ec21=A(upe)';
    it=0:-1:1-n;
    e12=rem(floor(ec12(:)*pow2(it)),2);
    e21=rem(floor(ec21(:)*pow2(it)),2);
    E=e21-e12;
    E(m,:)=ones(1,n);
    a=(v(ec21)-v(ec12))';
    a(m)=v(N);
    if n==2, a=a'; end;
    err=norm(E*x-a)^2; if err<eps, break; end
    Q=2*E'*E;
    b=sparse(2*E'*a);

% Calling conic programming solver.
    Q1=[E(m,:),zeros(1,n);Q,-eye(n)];
    Q1(:,end+1)=0;
    prob.a=sparse(Q1);
    prob.blc=[a(m);b];
    prob.buc=[a(m);b];
    lb=[-inf(2*n+1,1)];
    prob.blx=lb;
    c=[zeros(n,1);zeros(n,1);1];
    prob.c=c;
    prob.cones.type   = [res.symbcon.MSK_CT_QUAD];
    prob.cones.sub    = [2*n+1,n+1:2*n];
    prob.cones.subptr = [1];

    [r,res2] = mosekopt('minimize echo(0)',prob,param);
    sol=res2.sol;
    obj=sol.itr.pobjval;
    if strcmp(sol.itr.solsta,'OPTIMAL') ~= 1
       x=sol.itr.xx(1:n)';
       warning('ker:No','Probably no pre-kernel point found!');
       break;
    end
    x=sol.itr.xx(1:n);


% Due to a badly conditioned matrix, we might get an overflow/underflow.
% In this case, we restart with a new starting point.
    z1=any(isinf(x));
    z2=any(isnan(x));
    if z1==1 || z2==1 
       x=eye(n,1); 
    else 
    end
    Lerr(cnt,:)=[err, norm(E*x-a)^2]; % checking purpose
    xarr(cnt,:)=x'; % intermediate results
end

if cnt==CNT, % should trigger errors ....
  if slv==0 && smc==1
       msg01='No Pre-Kernel Element found. Changing Cardinality.';
       warning('PrK:ChangCard',msg01);
       if mnQ==1 && n < 15;x=4*x;end
       [x, Lerr, smat, xarr]=computePrk(v,x,0,slv,mnQ);
  elseif slv==0 && smc==0
       msg02='No Pre-Kernel Element found. Changing the Solver.';
       warning('PrK:ChangSolv',msg02);
       if mnQ==1; x=2*x'; else; x=LS_PreNucl(v); end
       [x, Lerr, smat, xarr]=computePrk(v,x,smc,1,mnQ);
  elseif slv==1 && smc==0
       msg01='No Pre-Kernel Element found. Changing Cardinality to Default Value.';
       warning('PrK:Default',msg01);
       if mnQ==1 && n < 7 || n>=11;x=(v(N)/n)*ones(n,1); end
       [x, Lerr, smat, xarr]=computePrk(v,x',1,1,mnQ);
  else
       x=x';
       msg02='No Pre-Kernel Element found. Change payoff vector and restart!';
       warning('PrK:NotFound',msg02);
  end
else
%  x=x';
end


%--------------
function [A, smat]=effCoalitions(v,x,smc,cnt)
% Computes the set of most effective coalitions
% of smallest/largest cardinality.
%
% Define variables:
% output:
% A     -- matrix of most effective coalitions of smallest/largest cardinality.
% smat  -- as above.
% cnt   -- loop counter.
%
% input:
% cnt   -- loop counter.
%       -- otherwise, as above.
%
n=length(x);
% The set of effective coalitions might be too
% large or too small due to floating point arithmetic.
% Adjusting the tolerance value might help to find the
% correct choice. In case that the set of most effective
% coalitions is not selected correctly, pathological
% cycles may appear.

if cnt<6
 tol=eps;
elseif cnt > 10
 tol=1500*eps;
else
 tol=100*eps;
end

% Borrowed from J. Derks
Xm=x(1); for ii=2:n, Xm=[Xm x(ii) Xm+x(ii)]; end
% Computing the excess vector w.r.t. x.
e=v-Xm;
clear v Xm;
% Truncate data arrays.
[e, sC]=sort(e,'descend');
B=eye(n);
smat=-inf(n);
q0=n^2-n;
q=0;
k=1;
pl=1:n;
while q~=q0
  kS=sC(k);
  ai=bitget(kS,pl)==1;
  bj=ai==0;
  pli=pl(ai);
  plj=pl(bj);
  if isempty(plj)==0
    for i=1:length(pli)
      for j=1:length(plj)
        if B(pli(i),plj(j))==0 
           B(pli(i),plj(j))=k;
           smat(pli(i),plj(j))=e(k); % max surplus of i against j.
           q=q+1;
        end
      end
    end
  end
  k=k+1;
end
m=max(B(:));
e1=e(m)-tol;
le=e>=e1;
tS=sC(le);
te=e(le);
clear e sC;
A=eye(n);
% Selecting the set of most effective coalitions 
% having smallest/largest cardinality.
% Assigning the set of selected coalitions to 
% matrix A
parfor i=1:n
   a=bitget(tS,i)==1;
   for j=1:n
    if i~=j
       b=bitget(tS,j)==0;
       lij=a & b;
       c_ij=tS(lij);
       ve=te;
       ex_ij=ve(lij);
       abest_ij=abs(smat(i,j)-ex_ij)<tol;
       slc_cij=c_ij(abest_ij);
       lC=length(slc_cij);
       if lC==1
          A(i,j)=slc_cij;
       else
          binCell_ij=SortSets(slc_cij,n,lC,smc);
          if smc==1
             A(i,j)=binCell_ij(1);  % Selecting smallest cardinality.
             elseif smc==0
             A(i,j)=binCell_ij(end); % Selecting largest cardinality.
          else
             A(i,j)=binCell_ij(end);   % Selecting largest cardinality.
          end
       end
    end
   end
end

%-------------------------------
function Seff=SortSets(effij,n,bd,smc)
% Sorting the set of most effective
% coalitions with respect to their
% cardinality. Ascent ordering.
% Smallest coalitions are coming first.
  Pm=zeros(bd,n);
  for k=1:n, Pm(:,k) = bitget(effij,k);end
  ov=ones(n,1);
  clsize=Pm*ov;
  if smc==1
     mcl=min(clsize);
  else
     mcl=max(clsize);
  end
  eqm=find(clsize==mcl);
  lc=length(eqm);
  if lc~=bd
     effij=effij(eqm);
     Pm=Pm(eqm,:);
     clsize=clsize(eqm);
  end
  pwcl=clsize.^3;
  J=1:n;
  J=J(ones(lc,1),:);
  M=Pm.*J;
  M=M.^(1/2);
  clix=M*ov;
  clnb=clix.*pwcl;
  [~, ix]=sort(clnb);
  Seff=effij(ix);

function [x, Lerr, smat, xarr]=p_msk_weightedPreKernel(clv,pS,x)
% P_MSK_WEIGHTEDPREKERNEL computes from (v,pS,x) a weighted pre-kernel element using qpas.mex
% and Matlab's PCT.
% Source: Meinhardt, 2010.
%
% MSK-Solver: http://www.mosek.com/  
%
% Usage: [x Lerr smat xarr] = p_msk_weightedPreKernel(clv,pS,x)
%
% Define variables:
%  output:
%  x        -- A weighted pre-kernel element (output)
%              (default per capita pre-kernel).
%  Lerr     -- List of computed function values of hx and h. 
%  smat     -- Matrix of maximum weighted surpluses.
%  xarr     -- History of computed solution at each iteration step.
%
%  input:
%  clv      -- TuGame class object.
%  pS       -- A vector of weights of length 2^n-1 (optional). 
%  x        -- payoff vector of size(1,n) (optional)


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/06/2014        0.5             hme
%                

if nargin<3
  v=clv.tuvalues;
  N=clv.tusize;
  n=clv.tuplayers;
  gt=clv.tutype;
  stx=clv.tustpt;
  mnQ=clv.tumnQ;
  if isempty(stx)
    Si=clv.tuSi;
    mv=clv.tumv;
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
       k=1:n;
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
    x=stx;
    smc=1;
  end
  if nargin < 2
    pS='';
  end
else
  v=clv.tuvalues;
  gt=clv.tutype;
  mnQ=clv.tumnQ;
  smc=1;
end

if islogical(v)
   v=double(v);
end

[x, Lerr, smat, xarr]=computePrk(v,x,pS,smc,0,mnQ);
smat=tril(smat,-1)+triu(smat,1);

% Main function to compute a
% pre-kernel element.
%-----------------------------
function [x, Lerr, smat, xarr]=computePrk(v,x,pS,smc,slv,mnQ)
% 
%  output:  -- as above.
%
%  input:
%   v       -- as above.
%   x       -- as above.
%  pS       -- A vector of weights of length 2^n-1.
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

ra = reasonable_outcome(v)';
k=1:n;
vi=v(bitset(0,k))';
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if isempty(pS)
  int=1-n:1:0;
  S=1:N;
  mat=(rem(floor(S(:)*pow2(int)),2)==1)';
  clS=ones(1,n)*mat;
  pS=1./clS; % weights of coalitions (per capita)
end



[rcode,res] = mosekopt('param echo(0)');
param=res.param;
param.MSK_IPAR_INTPNT_BASIS   = 'MSK_ON';
param.MSK_IPAR_OPTIMIZER = 1;  % Using interior-point optimizer.
param.MSK_DPAR_INTPNT_TOL_DFEAS = 1.0000e-12;
param.MSK_DPAR_INTPNT_TOL_PFEAS = 1.0000e-12;
param.MSK_DPAR_INTPNT_TOL_MU_RED = 1.0000e-14;
if n < 22
  param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
else
  param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-10;
end
%param.MSK_DPAR_INTPNT_TOL_REL_STEP = 0.25;
%param=[];
%param.MSK_DPAR_INTPNT_TOL_REL_GAP

% Cycling may occur, so that we need an artificial halt
while cnt<CNT  
    cnt=cnt+1;
    [A, smat]=effCoalitions(v,x,smc,cnt,pS);
    upe=tril(upe,-1);
    etr12=A';
    ec12=etr12(upe)';
    ec21=A(upe)';
    it=0:-1:1-n;
    e12=rem(floor(ec12(:)*pow2(it)),2);
    e21=rem(floor(ec21(:)*pow2(it)),2);
    pS12=pS(ec12)';
    pS21=pS(ec21)';
    mD12=diag(pS12);
    mD21=diag(pS21);    
    E=mD21*e21-mD12*e12;
    E(m,:)=ones(1,n);
    a=(pS(ec21).*v(ec21)-pS(ec12).*v(ec12))';
    a(m)=v(N);
    if n==2, a=a'; end;
    err=norm(E*x-a)^2; if err<eps, x=x';break; end
    Q=2*E'*E;
    b=sparse(-2*E'*a);
    [qi,qj,qs] = find(tril(Q));


% Calling quadratic programming solver.
    prob.qosubi=qi;
    prob.qosubj=qj;
    prob.qoval=qs;
    prob.c=b;
    subi=ones(n,1);
    subj=1:n;
    valij=E(m,:);
    prob.a = sparse(subi,subj,valij);
    prob.blc=a(m);
    prob.buc=a(m);
    prob.blx=[];
    prob.bux=sparse(ra);
    [r,res] = mosekopt('minimize echo(0)',prob,param);
    sol=res.sol;
    obj=sol.itr.pobjval;
    if strcmp(sol.itr.solsta,'OPTIMAL') ~= 1 
       x=sol.itr.xx';
       warning('ker:No','Probably no pre-kernel point found!');
       break;
    end
    x=sol.itr.xx;
%x'

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
       [x, Lerr, smat, xarr]=computePrk(v,x',pS,0,slv,mnQ);
  elseif slv==0 && smc==0
       msg02='No Pre-Kernel Element found. Changing the Solver.';
       warning('PrK:ChangSolv',msg02);
       if mnQ==1; x=2*x'; else; x=LS_PreNucl(v); end
       [x, Lerr, smat, xarr]=computePrk(v,x,pS,smc,1,mnQ);
  elseif slv==1 && smc==0
       msg01='No Pre-Kernel Element found. Changing Cardinality to Default Value.';
       warning('PrK:Default',msg01);
       if mnQ==1 && n < 7 || n>=11;x=(v(N)/n)*ones(n,1); end
       [x, Lerr, smat, xarr]=computePrk(v,x',pS,1,1,mnQ);
  else
       x=x';
       msg02='No Pre-Kernel Element found. Change payoff vector and restart!';
       warning('PrK:NotFound',msg02);
  end
else
%  x=x';
end


%--------------
function [A, smat]=effCoalitions(v,x,smc,cnt,clS)
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
N=2^n-1;
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
eN=e(N);
e=e.*clS;
e(N)=eN;
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
    for i=1:numel(pli)
      for j=1:numel(plj)
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
lcl=length(tS);
te=e(le);
clear e sC;


slcCell=cell(n);
A=eye(n);


% Selecting the set of most effective coalitions 
% having smallest/largest cardinality.

parfor i=1:n
   a=bitget(tS,i)==1;
   for j=1:n
     if i<j
       b=bitget(tS,j)==0;
       lij=a & b;
       c_ij=tS(lij);
       ex_ij=te(lij);
       abest_ij=abs(smat(i,j)-ex_ij)<tol;
       slcCell{i,j}=c_ij(abest_ij);
      elseif i>j
       b=bitget(tS,j)==0;
       lij=a & b;
       c_ij=tS(lij);
       ex_ij=te(lij);
       abest_ij=abs(smat(i,j)-ex_ij)<tol;
       slcCell{i,j}=c_ij(abest_ij);
    end
   end
end

% Assigning the set of selected coalitions to 
% matrix A.
parfor i=1:n
  for j=1:n
   if A(i,j)== 0
      lC=length(slcCell{i,j});
     if lC==1
        A(i,j)=slcCell{i,j}; 
     else
         binCell_ij=SortSets(slcCell{i,j},n,lC,smc);
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

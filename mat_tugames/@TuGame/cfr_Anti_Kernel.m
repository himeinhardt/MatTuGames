function [x, Lerr, smat, xarr]=cfr_Anti_Kernel(clv,F,x)
% CFR_ANTI_KERNEL computes from (v,x) a Kernel element with coalition formation restrictions using
% Matlab's Optimization Toolbox.
%
% Usage: [x Lerr smat xarr] = clv.cfr_Anti_Kernel(F,x)
%
% Define variables:
%  output:
%  x        -- Anti-Kernel element (output) w.r.t. coalition formation restrictions F.
%  Lerr     -- List of computed function values of hx and h. 
%  smat     -- Matrix of minimum surpluses.
%  xarr     -- History of computed solution at each iteration step.
%
%  input:
%  clv      -- TuGame class object.
%  F        -- For instance, a characterization set for the nucleolus.
%              F must contain the grand coalition N.
%  x        -- payoff vector of size(1,n) (optional)


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/22/2017        0.9             hme
%                

if nargin<2
  v=clv.tuvalues;
  N=clv.tusize;
  n=clv.tuplayers;
  gt=clv.tutype;
  stx=clv.tustpt;
  if isempty(stx)
    Si=clv.tuSi;
    mv=clv.tumv;
    x=(mv-v(Si))/2;
    sx=sum(x);
     if sx>0
       x=x*mv/sx;
      elseif all(abs(x-0)<10^3*eps)==1
       k=1:n;
       x=(mv-v(bitset(0,k)))/2;
       sx=sum(x);
       if sx>0
          x=select_starting_pt(clv);
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
else
  v=clv.tuvalues;
  gt=clv.tutype;
  smc=1;
end

if islogical(v)
   v=double(v);
end
si=clv.tuvi;
%% Addining the singleton coalitions to F, 
%% since v is definded over F ans si.
F=unique([F,si]);
[x, Lerr, smat, xarr]=cfr_computeAkr(v,x,smc,0,F);
smat=tril(smat,-1)+triu(smat,1);

% Main function to compute a
% pre-kernel element.
%-----------------------------
function [x, Lerr, smat, xarr]=cfr_computeAkr(v,x,smc,slv,F)
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

tol=10^8*eps;
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

ofval=inf;
ra = smallest_amount(v)';
k=1:n;
Nk=N-2.^(k-1);
vi=v(Nk)';
%vi=v(bitset(0,k))';
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if sum(vi)<v(N)
   error('sum of lower bound exceeds value of grand coalition! No solution can be found that satisfies the constraints.');
   return;
end

S=1:N;
lf=length(F);
lfNq=F(end)~=N;
if lfNq
   F(end+1)=N;
   lf=lf+1;
end
CS=S(ismember(S,F)==0);
vF=v;
vF(CS)=[];
if lfNq
   vF(end)=0;
end
clear S CS;

% Cycling may occur, so that we need an artificial halt
while cnt<CNT  
    cnt=cnt+1;
    [A, smat]=cfr_effCoalitions(v,x,smc,cnt,F);
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
    a(m)=vF(lf);
    if n==2, a=a'; end;
    err=norm(E*x-a)^2; if err<eps, x=x';break; end
% checking anti-kernel property
    ir=(x-vi)';
    irQ=all(ir>-tol);
    if irQ
      smat=tril(smat,-1)+triu(smat,1);
      krm=smat-smat';
      irm=repmat(ir,n,1);
      kriQ=all((krm.*irm)<=tol);
      effQ=abs(v(end)-sum(x))<tol;
      krQ=all(kriQ) && effQ;
    else
      krQ=false;
    end
    if krQ == 1; x=x'; break; end

    Q=2*E'*E;
    b=-2*E'*a;
    ra=[];
% Calling quadratic programming solver.
% Uncomment these lines if you don't have Mosek
    opts = optimset('Algorithm','interior-point-convex','Display','off','TolFun',1e-14);
% Comment the two lines about out if you don't have Mosek.
%    opts = optimset;
%    opts = optimset(opts,'MSK_DPAR_INTPNT_TOL_DFEAS',1.0000e-10,'MSK_DPAR_INTPNT_TOL_PFEAS',1.0000e-10,'MSK_DPAR_INTPNT_TOL_MU_RED',1.0000e-11,'MSK_DPAR_INTPNT_TOL_REL_GAP',1.0000e-11);
%    opts = optimset(opts,'MSK_DPAR_INTPNT_TOL_DFEAS',1.0000e-14,'MSK_DPAR_INTPNT_TOL_PFEAS',1.0000e-14,'MSK_DPAR_INTPNT_TOL_MU_RED',1.0000e-14,'MSK_DPAR_INTPNT_TOL_REL_GAP',1.0000e-14);

    [x,fval,exitflag,output,lambda] = quadprog(Q,b,[],[],E(m,:),a(m),ra,vi,x,opts);
    if exitflag ~= 1
       x=x';
       warning('Aker:No','Probably no anti-kernel point found!');
       break;
    elseif abs(fval-ofval)<tol
       if irQ
          smat=tril(smat,-1)+triu(smat,1);
          krm=smat-smat';
          irm=repmat(ir,n,1);
          kriQ=all((krm.*irm)<=tol);
          krQ=all(kriQ);
       else
          krQ=false;
       end
       x=x';
       if krQ==0
          warning('Aker:NoB','Probably no anti-kernel point found!');
       end
       break;
    end

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
       msg01='No Kernel Element found. Changing Cardinality.';
       warning('Kr:ChangCard',msg01);
       [x, Lerr, smat, xarr]=computeAkr(v,x',0,slv);
  else
       if irQ
          smat=tril(smat,-1)+triu(smat,1);
          krm=smat-smat';
          irm=repmat(ir,n,1);
          kriQ=all((krm.*irm)<=tol);
          krQ=all(kriQ);
       else
          krQ=false;
       end
       x=x';
       if krQ==0
          msg02='No Anti-Kernel Element found. Change payoff vector and restart!';
          warning('AKr:NotFound',msg02);
       end
  end
else
%  x=x';
end


%--------------
function [A, smat]=cfr_effCoalitions(v,x,smc,cnt,F)
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

%% F should contain the grand coalition for defining vF.
S=1:N;
CS=S(ismember(S,F)==0);
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
e(CS)=[];
clear v Xm;
% Truncate data arrays.
[e, sC]=sort(e,'ascend');
sC=F(sC);
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
           smat(pli(i),plj(j))=e(k); % min surplus of i against j.
           q=q+1;
        end
      end
    end
  end
  k=k+1;
end
m=max(B(:));
e1=e(m)+tol;
le=e<=e1;
tS=sC(le);
lcl=length(tS);
te=e(le);
clear e sC;

% Computing the set of most effective coalitions.
A=eye(n);
a=false(lcl,n);
c=cell(n);
slcCell=cell(n);
binCell=cell(n);
abest=cell(n);

% Constructing the set of coalitions containing player i
% without player j.
for i=1:n
   a(:,i)=bitget(tS,i)==1;
end
b=a==0;

% Selecting the set of most effective coalitions
% having smallest/largest cardinality.
for i=1:n-1
   for j=i+1:n
       lij=a(:,i) & b(:,j);
       lji=a(:,j) & b(:,i);
       c{i,j}=tS(lij);
       c{j,i}=tS(lji);
       ex_ij=te(lij);
       ex_ji=te(lji);
       abest{i,j}=abs(smat(i,j)-ex_ij)<tol;
       abest{j,i}=abs(smat(j,i)-ex_ji)<tol;
       slcCell{i,j}=c{i,j}(abest{i,j});
       slcCell{j,i}=c{j,i}(abest{j,i});
   end
end

% Assigning the set of selected coalitions to 
% matrix A.
for i=1:n-1
  for j=i+1:n
      lCi=length(slcCell{i,j});
      lCj=length(slcCell{j,i});
     if lCi==1
        A(i,j)=slcCell{i,j}; 
     else
         binCell{i,j}=SortSets(slcCell{i,j},n,lCi,smc);
      if smc==1
           A(i,j)=binCell{i,j}(1);  % Selecting smallest cardinality.
      elseif smc==0
           A(i,j)=binCell{i,j}(end); % Selecting largest cardinality.
      else
           A(i,j)=binCell{i,j}(end);   % Selecting largest cardinality.
      end
     end
     if lCj==1
        A(j,i)=slcCell{j,i};
     else
        binCell{j,i}=SortSets(slcCell{j,i},n,lCj,smc);
       if smc==1
           A(j,i)=binCell{j,i}(1);  % Selecting smallest cardinality.
       elseif smc==0
           A(j,i)=binCell{j,i}(end); % Selecting largest cardinality.
       else
           A(j,i)=binCell{j,i}(end);   % Selecting largest cardinality.
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

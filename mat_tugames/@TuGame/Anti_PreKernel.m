function [x, Lerr, smat, xarr]=Anti_PreKernel(clv,x)
%ANTI_PREKERNEL computes from (v,x) an anti-pre-kernel element.
% Source: Meinhardt, 2010.
%         Funaki and Meinhardt, 2006.
%
% Usage: [x Lerr smat xarr]=Anti_PreKernel(clv,x)
%
% Define variables:
%  output:
%  x        -- Anti-Pre-Kernel element (output)
%  Lerr     -- List of computed function values of hx and h.
%  smat     -- Matrix of minimum surpluses.
%  xarr     -- History of computed solution at each iteration step.
%
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n) (optional)

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/28/2012        0.3              hme
%   08/03/2016        0.9              hme
%   05/05/2019        1.1             hme
%

if nargin<2
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
      elseif all(abs(x)<10^3*eps)==1
       k=1:n;
       x=(mv-v(bitset(0,k)))/2;
       sx=sum(x);
       mmq=min(x)~=max(x);
       if sx>0 && mmq
          x=clv.select_starting_pt;
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
  N=clv.tusize;
  n=clv.tuplayers;
  gt=clv.tutype;
  mnQ=clv.tumnQ;
  smc=1;
end

gdata=struct('v',v,'x',x','smc',smc,'slv',0,'gt',gt,'mnQ',mnQ,'n',n,'N',N);

[x, Lerr, smat, xarr]=computeAntiPrk(gdata);
smat=tril(smat,-1)+triu(smat,1);

% Main function to compute a
% pre-kernel element.
%-----------------------------
function [x, Lerr, smat, xarr]=computeAntiPrk(gdata)
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

n=gdata.n;
N=gdata.N;
smc=gdata.smc;
mnQ=gdata.mnQ;
slv=gdata.slv;



cnt=0;
if 15<=n
 if mnQ==1
    CNT=n+9;
 else
    CNT=n+2;
 end
else
 CNT=2*(n+1);
end
Lerr=-inf(CNT,2);
xarr=-inf(CNT,n);
m=1+n*(n-1)/2;
upe=true(n);

% Cycling may occur, so that we need an artificial halt
while cnt<CNT  
    cnt=cnt+1;
    [A, smat]=effCoalitions(gdata,cnt);
    upe=tril(upe,-1);
    etr12=A';
    ec12=etr12(upe)';
    ec21=A(upe)';
    it=0:-1:1-n;
    e12=rem(floor(ec12(:)*pow2(it)),2);
    e21=rem(floor(ec21(:)*pow2(it)),2);
    E=e21-e12;
    E(m,:)=ones(1,n);
    a=(gdata.v(ec21)-gdata.v(ec12))';
    a(m)=gdata.v(N);
    if n==2, a=a'; end;
    err=norm(E*gdata.x-a)^2; if err<eps, break; end
    Q=E'*E;
    b=E'*a;
%
% If the set effc is unique and Q is non-singular all 
% solvers will compute the same solution.
  if n < 11 || islogical(gdata.v)
    if n < 7
        mtol=(5*n)^n*eps(max(diag(Q)));
        ctol=10^2;
    elseif n >= 7 && n <  9
        mtol=(5*n)^n*eps(max(diag(Q)));
        ctol=10^3;
    else
        mtol=30^n*eps(max(diag(Q)));
        ctol=10^3;
    end
    if slv==0
       if cond(Q)>1/mtol
           x1=qrginv(Q)*b; % (QR-decomposition)
           x2=svd_dec(Q,b); % (SVD-decomposition)
           hx1=norm(E*x1-a)^2;
           hx2=norm(E*x2-a)^2;
            if abs(hx2-hx1)<eps;
              gdata.x=x1;
            else
              gdata.x=x2;
            end
       else
% Q is symmetric and positive definite. We apply a Cholesky decomposition.
         gdata.x=chol_dec(Q,b); % (Cholesky factorization)
       end
    elseif slv==1
              if smc==0
                if cond(Q)>1/mtol
                 gdata.x=qr_dec(Q,b); % (QR-decomposition)
                else
                 if mnQ==0
                   gdata.x=svd_dec(Q,b); % (SVD-decomposition)
                 else
                   gdata.x=pinv(Q)*b;
                 end
                end
              elseif smc==1
                if cond(Q)>ctol
                 gdata.x=svd_dec(Q,b);  % (SVD-decomposition)
                else
                 gdata.x=chol_dec(Q,b); % (Cholesky factorization)
                end
              else
                if rcond(Q) < 1e-8
                 gdata.x=qr_dec(Q,b); % (QR-decomposition)
                else
                 gdata.x=qrginv(Q)*b;
                end
             end
    else 
        gdata.x=qrginv(Q)*b; % (QR-decomposition)
    end
  else
    if mnQ==1
       rcQ=rcond(Q);
       if rcQ > 1e-8
          gdata.x=svd_dec(Q,b);
       else
          if n<14
             Q=10^5*Q;
             b=10^5*b;
             gdata.x=qrginv(Q)*b;
          else
             gdata.x=pinv(E)*a;
          end
       end
    else
      gdata.x=qrginv(Q)*b;
    end
  end
%
% Due to a badly conditioned matrix, we might get an overflow/underflow.
% In this case, we restart with a new starting point.
    z1=any(isinf(gdata.x));
    z2=any(isnan(gdata.x));
    if z1==1 || z2==1 
       gdata.x=eye(n,1); 
    else 
    end
    Lerr(cnt,:)=[err, norm(E*gdata.x-a)^2]; % checking purpose
    xarr(cnt,:)=gdata.x'; % intermediate results
end

if err<eps
  x=gdata.x';
else
if cnt==CNT, % should trigger errors ....
  if slv==0 && smc==1
       msg01='No Anti Pre-Kernel Element found. Changing Cardinality.';
       warning('AnPrK:ChangCard',msg01);
       if mnQ==1 && n < 15;gdata.x=4*gdata.x;end 
       gdata.smc=0;
       [x, Lerr, smat, xarr]=computeAntiPrk(gdata);
  elseif slv==0 && smc==0 
       msg02='No Anti Pre-Kernel Element found. Changing the Solver.';
       warning('AnPrK:ChangSolv',msg02);
       if mnQ==1; gdata.x=2*gdata.x; else x=LS_PreNucl(gdata.v); end
       gdata.slv=1;  
       [x, Lerr, smat, xarr]=computeAntiPrk(gdata);
  elseif slv==1 && smc==0 
       msg01='No Anti Pre-Kernel Element found. Changing Cardinality to Default Value.';
       warning('AnPrK:Default',msg01);
       if mnQ==1 && n < 7 || n>=11;gdata.x=(gdata.v(N)/n)*ones(n,1); end 
       gdata.slv=1;
       gdata.smc=1;
       [x, Lerr, smat, xarr]=computeAntiPrk(gdata);
  elseif slv==1 && smc==1
       gdata.x=(gdata.v(N)/n)*ones(n,1);
       msg01='No Anti Pre-Kernel Element found. Changing to Start Value.';
       warning('AnPrK:StartVal',msg01);
       gdata.slv=1;
       gdata.smc=2;
       [x, Lerr, smat, xarr]=computeAntiPrk(gdata);
  else
       x=gdata.x';
       msg02='No Anti Pre-Kernel Element found. Change payoff vector and restart!';
       warning('AnPrK:NotFound',msg02);
  end
else
  x=gdata.x';
end
end
%--------------------------
function x=chol_dec(Q,b)
% Cholesky factorization solves the system
% Q x = b whenever Q is symmetric and positive definite.
%
if rcond(Q) < 5e-05
   Q=sparse(Q);
end
R=chol(Q);
y=R'\b;
x=R\y;
%-----------------------------
function x=qr_dec(Q,b)
% QR-decomposition in order to solve the system 
% Q x = b 
% 
%
mf=factorize(Q);
if mf.A_condest > 10^16
  if rcond(Q) < 5e-5
     Q=sparse(Q);
     [Q1,R1]=spqr(Q);
  else
     [Q1,R1]=qr(Q,0);
  end
  x=Q1*(R1'\b);
else
   x=mf\b;
end
%------------------------
function x=svd_dec(Q,b)
% SVD-decomposition in order to solve the
% system Q x = b.
%
n=size(Q);
if n < 9
  th=1e-18;
else
  th=1e-17;
end
if rcond(Q) < th
  F=linfactor(Q);
  x=linfactor(F,b);
else
  [U1,S1,V]=svd(Q,0);
  y=S1\(U1'*b);
  x=V*y;
end

%--------------
function [A, smat]=effCoalitions(gdata,cnt)
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
% gt    -- game type.
%       -- otherwise, as above.
%
n=gdata.n;
smc=gdata.smc;
mnQ=gdata.mnQ;
N=gdata.N;
% The set of effective coalitions might be too
% large or too small due to floating point arithmetic.
% Adjusting the tolerance value might help to find the
% correct choice. In case that the set of most effective
% coalitions is not selected correctly, pathological
% cycles may appear.

if cnt<6
 if strcmp(gdata.gt,'cv')
  tol=5000*eps;
 else
  tol=eps;
 end
elseif cnt > 10
 tol=1500*eps;
else
 tol=100*eps;
end
%
% Inspired by Jean Derks.
Xm{1}=gdata.x(1); for ii=2:n, Xm{1}=[Xm{1} gdata.x(ii) Xm{1}+gdata.x(ii)]; end
% Computing the excess vector w.r.t. x.
e=gdata.v-Xm{1};
[e,sC]=sort(e,'ascend');
% Truncate data arrays.
B=eye(n);
smat=-inf(n);
q0=n^2-n;
q=0;
k=1;
pl=1:n;
while q~=q0
  kS=sC(k);
  if kS<N
    ai=bitget(kS,pl)==1;
    bj=ai==0;
    pli=pl(ai);
    plj=pl(bj);
    for i=1:length(pli)
      for j=1:length(plj)
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
te=e(le);
clear e sC;
% Computing the set of most effective coalitions.
A=eye(n);
% Constructing the set of coalitions containing player i
% without player j.
% Selecting the set of less effective coalitions
% having smallest/largest cardinality.
for i=1:n
   a=bitget(tS,i)==1;
   for j=1:n
    if i~=j
       b=bitget(tS,j)==0;
       lij=a & b;
       cij=tS(lij);
       ex_ij=te(lij);
       abest=abs(smat(i,j)-ex_ij)<tol;
       slc_cij=cij(abest);
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

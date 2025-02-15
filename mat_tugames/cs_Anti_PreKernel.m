function [x, Lerr, smat, xarr]=cs_Anti_PreKernel(v,cs,x)
%CS_ANTI_PREKERNEL computes from (v,x) an anti-pre-kernel element w.r.t. the coalition structure cs.
%
% Usage: [x Lerr smat xarr]=cs Anti_PreKernel(v,cs,x)
%
% Define variables:
%  output:
%  x        -- Anti-Pre-Kernel element (output) w.r.t. the coalition structure cs.
%  Lerr     -- List of computed function values of hx and h.
%  smat     -- Matrix of minimum surpluses.
%  xarr     -- History of computed solution at each iteration step.
%  input:
%  v       -- The dual game of v of length 2^n-1.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  x        -- payoff vector of size(1,n) (optional)

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/02/2016        0.9              hme
%   05/15/2022        1.9.1           hme
%


if nargin<1
    error('At least the game must be given!');
elseif nargin<2 
    error('A coalition structure cs must be provided!');
elseif nargin<3
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
     elseif all(abs(x)<10^3*eps)==1
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
x=x';
cs=double(cs);
[x, Lerr, smat, xarr]=cs_computeAntiPrk(v,x,smc,0,mnQ,cs);
smat=tril(smat,-1)+triu(smat,1);

% Main function to compute a
% pre-kernel element.
%-----------------------------
function [x, Lerr, smat, xarr]=cs_computeAntiPrk(v,x,smc,slv,mnQ,cs)
% 
%  output:  -- as above.
%
%  input:
%   v       -- as above.
%   cs      -- as above.
%   x       -- as above.
%  smc      -- selecting from effc the smallest/largest cardinality (optional).
%              Value must be set to 0 (largest),1 (smallest) or 2 (reset).
%  slv      -- selecting a different linear solver (QR/SVD-decomposition). 
%              Value must be set to 0 or 1.

n=length(x);
N=2^n-1;
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
lcs=length(cs);
mb=1:n;
for kk=1:lcs 
    clm{kk}=mb(bitget(cs(kk),mb)==1);
end

% Cycling may occur, so that we need an artificial halt
while cnt<CNT  
    cnt=cnt+1;
    [A, smat]=cs_effCoalitions(v,x,smc,cnt);
    if length(cs)>1
       B12=-inf(n,n);
       B21=-inf(n,n);
       for kk=1:lcs
           ckk=clm{kk};
           lckk=numel(clm{kk});
           if lckk>1 
              for ii=1:lckk-1
                   for jj=ii+1:lckk
                       B12(ckk(ii),ckk(jj))=A(ckk(ii),ckk(jj));
                       B21(ckk(jj),ckk(ii))=A(ckk(jj),ckk(ii));
                   end
              end
           end
       end
       bm12=B12(:);
       B21=B21';
       bm21=B21(:);
       gif12=bm12>-inf;
       gif21=bm21>-inf;
       bm12=bm12(gif12);
       bm21=bm21(gif21);
       it=0:-1:1-n;
       e12=rem(floor(bm12(:)*pow2(it)),2);
       e21=rem(floor(bm21(:)*pow2(it)),2);
       E=e21-e12;
       a=(v(bm21)-v(bm12))';
       acs=rem(floor(cs(:)*pow2(it)),2);
       E=[E;acs]; 
       [s1,s2]=size(E);
       a(end+1:s1,1)=v(cs)';
    else   
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
    end

    if n==2, a=a'; end;
    err=norm(E*x-a)^2; 
    if err<eps, break; end
    Q=E'*E;
    b=E'*a;
%
% If the set effc is unique and Q is non-singular all 
% solvers will compute the same solution.
  if n < 11 || islogical(v)
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
              x=x1;
            else
              x=x2;
            end
       else   
% Q is symmetric and positive definite. We apply a Cholesky decomposition.
	 x=chol_dec(Q,b); % (Cholesky factorization)
       end
    elseif slv==1
              if smc==0 
                if cond(Q)>1/mtol
                 x=qr_dec(Q,b); % (QR-decomposition)
                else
                 if mnQ==0
	           x=svd_dec(Q,b); % (SVD-decomposition)
                 else
                   x=pinv(Q)*b;
                 end
                end
              elseif smc==1 
                if cond(Q)>ctol
                   x=svd_dec(Q,b);  % (SVD-decomposition)
                else 
	           x=chol_dec(Q,b); % (Cholesky factorization)
                end
              else 
                if rcond(Q) < 1e-8
                 x=qr_dec(Q,b); % (QR-decomposition)
                else
                 x=qrginv(Q)*b;
                end   
             end
    else 
        x=qrginv(Q)*b; % (QR-decomposition)
    end
  else
    if mnQ==1
       rcQ=rcond(Q);
       if rcQ > 1e-8
          x=svd_dec(Q,b);
       else
          if n<14
             Q=10^5*Q;
             b=10^5*b;
             x=qrginv(Q)*b;
          else
             x=pinv(E)*a;
          end
       end
    else
     x=qrginv(Q)*b;
    end
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

if err<eps
  x=x';
else
if cnt==CNT, % should trigger errors ....
  if slv==0 && smc==1
       msg01='No Anti-Pre-Kernel Element found. Changing Cardinality.';
       warning('AnPrK:ChangCard',msg01);
       [x, Lerr, smat, xarr]=cs_computeAntiPrk(v,x,0,slv,mnQ,cs);
  elseif slv==0 && smc==0 
       msg02='No Anti-Pre-Kernel Element found. Changing the Solver.';
       warning('AnPrK:ChangSolv',msg02);
       [x, Lerr, smat, xarr]=cs_computeAntiPrk(v,x,smc,1,mnQ,cs);
  elseif slv==1 && smc==0 
       msg01='No Anti-Pre-Kernel Element found. Changing Cardinality to Default Value.';
       warning('AnPrK:Default',msg01);
       [x, Lerr, smat, xarr]=cs_computeAntiPrk(v,x,1,1,mnQ,cs);
  elseif slv==1 && smc==1
       x=(v(N)/n)*ones(1,n);
       msg01='No Anti-Pre-Kernel Element found. Changing to Start Value.';
       warning('AnPrK:StartVal',msg01);
       [x, Lerr, smat, xarr]=cs_computeAntiPrk(v,x',2,1,mnQ,cs);
  else
       x=x';
       msg02='No Anit-Pre-Kernel Element found. Change payoff vector and restart!';
       warning('AnPrK:NotFound',msg02);
  end
else
  x=x';
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
%
%--------------
function [A, smat]=cs_effCoalitions(v,x,smc,cnt)
% Computes the set of most effective coalitions
% of smallest/largest cardinality.
%
% Define variables:
% output:
% A     -- matrix of most effective coalitions of smallest/largest cardinality.
% smat  -- as above.
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

% Inspired by Jean Derks.
Xm{1}=x(1); for ii=2:n, Xm{1}=[Xm{1} x(ii) Xm{1}+x(ii)]; end
% Computing the excess vector w.r.t. x.
e=v-Xm{1};
clear v Xm;
% Truncate data arrays.
[e, sC]=sort(e,'ascend');
B=zeros(n,n);
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
       abest_ij=abs(smat(i,j)-ex_ij)<tol;
       abest_ji=abs(smat(j,i)-ex_ji)<tol;
       slcCell{i,j}=c{i,j}(abest_ij);
       slcCell{j,i}=c{j,i}(abest_ji);
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
% Those of equal size are order lexicographically.     
pl=1:n;
bd=length(effij);
it=0:-1:1-n;
indM=rem(floor(effij(:)*pow2(it)),2);
ov=ones(n,1);
clsize=indM*ov;
  if smc==1
     mcl=min(clsize);
  else
     mcl=max(clsize);
  end
  eqm=find(clsize==mcl);
  lc=length(eqm);
  if lc~=bd
     effij=effij(eqm);
     indM=indM(eqm,:);
  end
expl=pl.*ones(lc,n);
imat=indM.*expl;
clm=n-imat;
clm=2.^clm;
dln=clm*ov;
%dln=sum(clm,2)'; %% canonical numbers of coalitions S.
%% canonical order R lex T iff d(R) > d(T).
[st,sid]=sort(dln,'descend');
%% Determining canonical order of coalitions S.
Seff=effij(sid);

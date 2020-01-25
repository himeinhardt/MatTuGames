function [x, Lerr, smat, xarr]=Anti_PreKernel(v,x)
%ANTI_PREKERNEL computes from (v,x) an anti-pre-kernel element.
% Source: Meinhardt, 2010.
%         Funaki and Meinhardt, 2006.
%
% Usage: [x Lerr smat xarr]=Anti_PreKernel(v,x)
%
% Define variables:
%  output:
%  x        -- Anti-Pre-Kernel element (output)
%  Lerr     -- List of computed function values of hx and h.
%  smat     -- Matrix of minimum surpluses.
%  xarr     -- History of computed solution at each iteration step.
%  input:
%  v       -- The dual game of v of length 2^n-1.
%  x        -- payoff vector of size(1,n) (optional)

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/16/2012        0.2 beta 1       hme
%   09/12/2012        0.2              hme
%   10/25/2012        0.3              hme
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
    Si=bitset(N,k,0);
    mv=max(v);
    x=(mv-v(Si))/2;
    sx=sum(x);
     if sx>0
       x=x*mv/sx;
      elseif all(abs(x-0)<10^3*eps)==1
       x=(mv-v(bitset(0,k)))/2;
       sx=sum(x);
       if sx>0
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
    if (2^n-1)~=N
       error('Game has not the correct size!'); 
    end
    smc=1;
end

[x, Lerr, smat, xarr]=computeAntiPrk(v,x,smc,0);
smat=tril(smat,-1)+triu(smat,1);

% Main function to compute a
% pre-kernel element.
%-----------------------------
function [x, Lerr, smat, xarr]=computeAntiPrk(v,x,smc,slv)
% 
%  output:  -- as above.
%
%  input:
%   v       -- as above.
%   x       -- as above.
%  smc      -- selecting from effc the smallest/largest cardinality (optional).
%              Value must be set to 0 (largest),1 (smallest) or 2 (reset).
%  slv      -- selecting a different linear solver (QR/SVD-decomposition). 
%               Value must be set to 0 or 1.

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
    Q=E'*E;
    b=E'*a;
%
% If the set effc is unique and Q is non-singular all 
% solvers will compute the same solution.
    if slv==0
       if cond(Q)>10^4
           x1=pinv(Q)*b; %x=pinv(E)*a (SVD-decomposition)
           x2=svd_dec(Q,b);
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
                if cond(Q)>10^3
                 x=qr_dec(Q,b); % (QR-decomposition)
                else 
	         x=chol_dec(Q,b); % (Cholesky factorization)               
                end
              elseif smc==1 
                if cond(Q)>10^3
                 x=svd_dec(Q,b);  % (SVD-decomposition)
                else 
	         x=chol_dec(Q,b); % (Cholesky factorization)               
                end
              else 
                if cond(Q)>10^3
                 x=qr_dec(Q,b); % (QR-decomposition)
                else 
	         x=chol_dec(Q,b); % (Cholesky factorization)               
                end   
             end
    else x=pinv(Q)*b; % (SVD-decomposition)
    end
    if cnt>20
        cnt;
    end
% Due to a badly conditioned matrix, we might get an overflow/underflow.
% In this case, we restart with a new starting point.
    z1=any(isinf(x));
    z2=any(isnan(x));
    if z1==1 || z2==1 
       idm=eye(n); 
       x=idm(:,1); 
    else 
    end
    Lerr(cnt,:)=[err, norm(E*x-a)^2]; % checking purpose
    xarr(cnt,:)=x'; % intermediate results
end

if cnt==CNT, % should trigger errors ....
  if slv==0 && smc==1
       msg01='No Anti-Pre-Kernel Element found. Changing Cardinality.';
       warning('AnPrK:ChangCard',msg01);
       [x, Lerr, smat, xarr]=computeAntiPrk(v,x',0,slv);
  elseif slv==0 && smc==0 
       msg02='No Anti-Pre-Kernel Element found. Changing the Solver.';
       warning('AnPrK:ChangSolv',msg02);
       [x, Lerr, smat, xarr]=computeAntiPrk(v,x',smc,1);
  elseif slv==1 && smc==0 
       msg01='No Anti-Pre-Kernel Element found. Changing Cardinality to Default Value.';
       warning('AnPrK:Default',msg01);
       [x, Lerr, smat, xarr]=computeAntiPrk(v,x',1,1);
  elseif slv==1 && smc==1
       x=(v(N)/n)*ones(1,n);
       msg01='No Anti-Pre-Kernel Element found. Changing to Start Value.';
       warning('AnPrK:StartVal',msg01);
       [x, Lerr, smat, xarr]=computeAntiPrk(v,x,2,1);
  else
       x=x';
       msg02='No Anit-Pre-Kernel Element found. Change payoff vector and restart!';
       warning('AnPrK:NotFound',msg02);
  end
else
  x=x';
end

%--------------------------
function x=chol_dec(Q,b);
% Cholesky factorization solves the system
% Q x = b whenever Q is symmetric and positive definite.
%
R=chol(Q);
y=R'\b;
x=R\y;

%-----------------------------
function x=qr_dec(Q,b);
% QR-decomposition in order to solve the system 
% Q x = b 
% 
%
[Q1,R1,P1]=qr(Q);
y=R1'\b;
x=Q1*y;


%------------------------
function x=svd_dec(Q,b);
% SVD-decomposition in order to solve the
% system Q x = b.
%
[U1,S1,V]=svd(Q);
y=S1\(U1'*b);
x=V*y;

%--------------
function [A, smat]=effCoalitions(v,x,smc,cnt);
% Computes the set of most effective coalitions
% of smallest/largest cardinality.
%
% Define variables:
% output:
% A     -- matrix of most effective coalitions of smallest/largest cardinality.
% smat  -- as above.
%
% input: 
%       -- as above.
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
v=[];
Xm=[];
% Truncate data arrays. 
[e, sC]=sort(e,'ascend');
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
e=[];
sC=[];

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
       c{i,j}=tS(a(:,i) & b(:,j));
       c{j,i}=tS(a(:,j) & b(:,i));
       ex_ij=te(a(:,i) & b(:,j));
       ex_ji=te(a(:,j) & b(:,i));
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
function Seff=SortSets(effij,n,bd,smc);
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

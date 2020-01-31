function [x, Lerr, smat, xarr]=CddPreKernel(v,x)
% CDDPREKERNEL computes from (v,x) a pre-kernel element.
% Inspired by Jean Derks, see email from 26/05/2010.
% Source: Meinhardt, 2010.
% 
%
% Usage: [x Lerr smat xarr]=CddPreKernel(v,x)
% Define variables:
%  output:
%  x        -- Pre-Kernel element (output)
%  Lerr     -- List of computed function values of hx and h. 
%  smat     -- Matrix of maximal surpluses.
%  xarr     -- History of computed solution at each iteration step.
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
%   07/10/2010        0.1 alpha       hme
%   02/13/2011        0.1 beta 2      hme
%   05/20/2012        0.2 beta        hme
%   10/25/2012        0.3 beta        hme
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
else
    N=length(v); 
    [~, n]=log2(N); 
    if (2^n-1)~=N
       error('Game has not the correct size!'); 
    end
end

CNT=2*(n+1);
x=x';
Lerr=[];
xarr=[];;
for i=1:CNT
  [x, Lerr01, smat, xarr01]=SeekPrk(v,x,1,n,N);
  Lerr=[Lerr; Lerr01];
  xarr=[xarr; xarr01];
end

smat=tril(smat,-1)+triu(smat,1);
if all(all(abs(smat-smat')<2*10^6*eps))==1,
   x=x';
else
 msg01='No Pre-Kernel Element found. Changing Cardinality.';
 warning(msg01);
 for i=1:CNT
   [x, Lerr01, smat, xarr01]=SeekPrk(v,x,0,n,N);
   Lerr=[Lerr; Lerr01];
   xarr=[xarr; xarr01];
 end
x=x';
smat=tril(smat,-1)+triu(smat,1);
msg02='No Pre-Kernel Element found. Change payoff vector and restart!';
warning(msg02);
end



%-----------------------------------------
function [x, Lerr, smat, xarr]=SeekPrk(v,x,smc,n,N)
%
%

% Define variables:
%  output:
%  x        -- solution of the linear equation Q*x = b.
%  Q        -- a symmetric and (semi)-positive definite matrix.
%  b        -- column vector b.
%  e        -- excess vector of game v w.r.t. x (input).
%  smat     -- matrix of maximal surpluses.
%  Ac       -- matrix of smallest/largest effective coalitions.
%  E        -- matrix of best arguments.
%  binCell  -- Cell contains the set of most effective coalitions.
%              binary representation.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)
%  smc      -- The cardinality of the set of most effective coalitions.
%              Smallest and largest cardinality can be chosen.
%              Permissible values are:
%              1 to invoke the samllest cardinality. 
%              0 to invoke the largest cardinality.
%

Lerr=[];
xarr=[];
m=1+n*(n-1)/2;
upe=true(n);


[Ac, smat]=effCoalitions(v,x,smc);
upe=triu(upe,1);
ec12=Ac(upe)';
etr21=Ac';
ec21=(etr21(upe))';
it=0:-1:1-n;
e12=rem(floor(ec12(:)*pow2(it)),2);
e21=rem(floor(ec21(:)*pow2(it)),2);
E=e21-e12;
E(m,:)=ones(1,n);
a=(v(ec21)-v(ec12))';
a(m)=v(N);
  if n==2, a=a'; end;
  err=norm(E*x-a)^2; if err<eps, return; end
Q=E'*E;
b=E'*a;
A1=Q;
B1=b;
H=struct('A',A1,'B',B1,'lin',(1:size(B1,1))');
V=cddmex('extreme',H);
x=V.V';

[rw cl]=size(x);
if cl>1
  y=pinv(Q)*b;
  x(:,cl+1)=y;
  for k=1:cl+1 nmi(k)=norm(x(:,k)); end
  [mv, id] = max(nmi);
  x=x(:,id);
else 
end
%x=pinv(Q)*b;
% Due to a badly conditioned matrix, we might get an overflow/underflow.
% In this case, we restart with a new starting point.
z1=any(isinf(x));
z2=any(isnan(x));
 if z1==1 || z2==1 
     idm=eye(n); 
     x=idm(:,1); 
   else 
  end
Lerr=[Lerr; err, norm(E*x-a)^2]; % checking purpose
xarr=[xarr; x']; % intermediate results


%-----------------------------
function x=qr_dec(Q,b,n)
% QR-decomposition in order to solve the system 
% Q x = b 
% 
%
[Q1, R1]=qr(Q,0);
y=R1'\b;
x=Q1*y;


%--------------
function [A, smat]=effCoalitions(v,x,smc)
% Computes the set of most effective coalitions
% of smallest/largest cardinality.
%
% Define variables:
% output:
% A       -- matrix of most effective coalitions of smallest/largest cardinality.
% smat    -- as above.
%
% input:  -- as above.


n=length(x);
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
tol=5000*eps;
e1=e(m)-tol;
le=e>=e1;
tS=sC(le);
lcl=length(tS);
te=e(le);
clear e sC;
% Computing the set of most effective coalitions.
A=eye(n);
% Constructing the set of coalitions containing player i
% without player j.
% Selecting the set of most effective coalitions
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
       lC=numel(slc_cij);
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

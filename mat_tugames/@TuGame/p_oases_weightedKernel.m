function [x, Lerr, smat, xarr]=p_oases_weightedKernel(clv,pS,x)
% P_OASES_WEIGHTEDKERNEL computes from (v,pS,x) a weighted kernel element using qpOASESmex
% and Matlab's PCT.
% Source: Meinhardt, 2010.
%
%  http://www.kuleuven.be/optec/software/qpOASES
%
%
% Usage: [x Lerr smat xarr] = p_oases_weightedKernel(clv,pS,x)
%
% Define variables:
%  output:
%  x        -- A weighted kernel element (output)
%              (default per capita kernel).
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
  essQ=clv.tuessQ;
  if essQ==0
     error('Game is not essential!')
  end
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
  essQ=clv.tuessQ;
  if essQ==0
     error('Game is not essential!')
  end
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
if isempty(pS)
  int=1-n:1:0;
  S=1:N;
  mat=(rem(floor(S(:)*pow2(int)),2)==1)';
  clS=ones(1,n)*mat;
  pS=1./clS; % weights of coalitions (per capita)
end

ofval=inf;
ra = reasonable_outcome(v)';
k=1:n;
vi=v(bitset(0,k))';
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end

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
% checking kernel property
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
      krQ=0;
    end
    if krQ == 1; x=x'; break; end

    Q=2*E'*E;
    b=-2*E'*a;
    A=E(m,:);
    lbA=a(m);
    ubA=lbA;

% Observing incorrect computation!
% Calling quadratic programming solver.
    opts = qpOASES_options;
%    opts.enableRegularisation=true;
    opts.enableCholeskyRefactorisation=true;
%    opts.enableEqualities=true;
    opts.boundTolerance=2.2204e-12;
    opts.terminationTolerance=2.2204e-12;
    [x,fval,exitflag,output,lambda] = qpOASES(Q,b,A,vi,ra,lbA,ubA,x,opts);
    if exitflag ~= 0
       x=x';
       warning('ker:No','Probably no kernel point found!');
       break;
    elseif abs(fval-ofval)<tol
       if irQ
          smat=tril(smat,-1)+triu(smat,1);
          krm=smat-smat';
          irm=repmat(ir,n,1);
          kriQ=all((krm.*irm)<=tol);
          krQ=all(kriQ);
       else
          krQ=0;
       end
       x=x';
       if krQ==0
          warning('Ker:NoB','Probably no kernel point found!');
       end
       break;
    end
    ofval=fval;

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
       if mnQ==1 && n < 15;x=4*x;end
       [x, Lerr, smat, xarr]=computePrk(v,x',pS,0,slv,mnQ);
  else
       if irQ
          smat=tril(smat,-1)+triu(smat,1);
          krm=smat-smat';
          irm=repmat(ir,n,1);
          kriQ=all((krm.*irm)<=tol);
          krQ=all(kriQ);
       else
          krQ=0;
       end
       x=x';
       if krQ==0
          msg02='No Kernel Element found. Change payoff vector and restart!';
          warning('Kr:NotFound',msg02);
       end
  end
else
%  warning('KerB:No','Probably no kernel point found!');
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

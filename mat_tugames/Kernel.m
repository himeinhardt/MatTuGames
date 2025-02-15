function [x, Lerr, smat, xarr]=Kernel(v,x)
% KERNEL computes from (v,x) a Kernel element using
% Matlab's Optimization Toolbox.
% Source: Meinhardt, 2013.
%
% Usage: [x Lerr smat xarr] = Kernel(v,x)
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
%   12/28/2012        0.3             hme
%   05/11/2019        1.1             hme
%   10/21/2020        1.9             hme
%   09/13/2021        1.9.1           hme
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
    vi=v(bitset(0,k));
    slb=sum(vi)>v(N);
    if slb==1
      error('Game is not essential!')
    end
    if N==1
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
    k=1:n;
    vi=v(bitset(0,k));
    slb=sum(vi)>v(N);
    if slb==1
      error('Game is not essential!')
    end 
    mv=max(v);
    mnQ=mv>v(N);
    if (2^n-1)~=N
       error('Game has not the correct size!'); 
    end
    smc=1;
end


[x, Lerr, smat, xarr]=computePrk(v,x,smc,0,mnQ);
smat=tril(smat,-1)+triu(smat,1);
end
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
ra=[];
k=1:n;
vi=double(v(bitset(0,k)))';


% Cycling may occur, so that we need an artificial halt
while cnt<CNT  
    cnt=cnt+1;
    [A, smat, pv]=effCoalitions(v,x,smc,cnt,[]);
    if isempty(A) %% is empty if truncating of the data array has failed.
       [A, smat]=effCoalitions(v,x,smc,cnt,pv);
    end
    if isempty(A)
       msg00='Probably No Pre-Kernel Element found!.';
       warning('PrK:ChangCard',msg00);
       x=-inf(1,n);
       Lerr=[];
       smat=[];
       xarr=[];
       return;
    end
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
    if n==2, a=a'; end
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

% Calling quadratic programming solver.
% Setting options
%    opts = optimset('Algorithm','interior-point-convex','Display','off','TolFun',1e-12);
    opts = optimoptions('quadprog','Algorithm','active-set','Display','off','TolFun',1e-12);
    ws=optimwarmstart(x',opts);
    [ws,fval,exitflag,output,lambda] = quadprog(Q,b,[],[],E(m,:),a(m),vi,ra,ws);
    if exitflag ~= 1
       x=ws.X';
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
       x=ws.X'
       if krQ==0
          warning('Ker:NoB','Probably no kernel point found!');
       end
       break;
    end
    x=ws.X;
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

if cnt==CNT % should trigger errors ....
  if slv==0 && smc==1
       msg01='No Kernel Element found. Changing Cardinality.';
       warning('Kr:ChangCard',msg01);
       if mnQ==1 && n < 15;x=4*x;end
       [x, Lerr, smat, xarr]=computePrk(v,x',0,slv,mnQ);
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
%  x=x';
end
end

%--------------
function [A, smat, pv]=effCoalitions(v,x,smc,cnt,pv)
% Computes the set of most effective coalitions
% of smallest/largest cardinality.
%
% Define variables:
% output:
% A     -- matrix of most effective coalitions of smallest/largest cardinality.
% smat  -- as above.
% pv    -- lower bound for the excesses that be treated.
%
% input:
% cnt   -- loop counter.
%       -- otherwise, as above.
%
n=length(x);
N=2^n-1;
if isempty(pv)
   pv1=0;
else
   pv1=pv;
end
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
lv=islogical(v);
%e=v-Xm{1};
e=bsxfun(@minus,v,Xm{1});
nf=max(randi([1 N],1,5));
el=e(nf);
% Truncate excess vector.
if n>16
   %el = min(e);
   eh = max(e);
   if abs(eh+el)<10^7*eps
      clear v Xm;
      [e,sC]=sort(e,'descend');
   elseif eh==1 && el==0
      clear v Xm;
      [e,sC]=sort(e,'descend');
   else
      if lv==1
         clear v Xm;
         if eh > 0.7 && eh < 1
            eh=1.3*eh;
            pv=min((eh+el)*0.9,0.7); % 0.6 fine
         else
            pv=min((eh+el)*0.3,0.8);
         end
      else
         vN=v(N);
         [mv,idx]=max(v);
         clear v Xm;
         if mv<0
            eh= 0.3*eh;
            pv=(eh+el)*1.2; % eh 0.3 fine (increase to improve).     
         elseif mv>vN
            eh = 0.8*eh;
            pv=(eh+el)*1.2; % 0.9 fine (increase to improve).
         elseif idx<N
            eh = 0.8*eh;
            pv=(eh+el); % 0.9 fine (increase to improve).
         else
            eh = 0.8*eh;
            pv=(eh+el)*0.3;
         end
      end
      if pv1==0
         pv=0.9*min(eh,pv);
      else
         pv=min(pv,pv1/2);
      end
      lp=e>pv-tol;
      e=e(lp);
      fS=find(lp);
      [e,fC]=sort(e,'descend');
      sC=fS(fC);
   end
else
   [e,sC]=sort(e,'descend');
end
% Truncate data arrays.
lC=length(sC);
B=eye(n);
smat=-inf(n);
q0=n^2-n;
q=0;
k=1;
pl=1:n;
while q~=q0
  kS=sC(k);
  if kS<N
     ai=(bsxfun(@bitget,kS,pl))==1;
     bj=ai==0;
     pli=pl(ai);
     for i=pli
        ri=B(i,:)==0;
        if any(ri)
           bn=pl(ri);
           plj=pl(bj);
           sj=plj(ismembc(plj,bn));
           lj=length(sj);
           if lj>0
             B(i,sj)=k;
             smat(i,sj)=e(k); % max surplus of i against j.
             q=q+lj;
           end
        end
     end
  end
  k=k+1;
  if k>lC
     break;
  end
end
%%% If we did too much concerning truncating
%%% then we restart from scratch.
if k>lC
   A=[];
   smat=[];
   return;
end
m=max(B(:));
e1=e(m)-tol;
le=e>=e1;
tS=sC(le);
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
end

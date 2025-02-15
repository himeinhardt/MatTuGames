function [krQ, kriQ, smat]=weightedCsKernelQ(clv,cs,x,pS,tol)
% WEIGHTEDCSKERNELQ checks whether the imputation x is a weighted kernel element 
% of the TU-game v w.r.t. coalition structure cs..
% 
%  Usage:[krQ, kriQ, smat]=clv.weightedCsKernelQ(cs,x,pS,tol);
%
%
% Define variables:
%  output:
%  krQ      -- Returns 1 (true) whenever the impuatation x is 
%              a weighted kernel element, otherwise 0 (false).
%  kriQ     -- Returns a list of true (1) and/or false (0) of length n.
%  smat     -- Matrix of maximum weighted surpluses.
%  input:
%  clv      -- TuGame class object.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  x        -- payoff vector of size(1,n) (optional)
%  pS       -- A vector of weights of length 2^n-1. (optional)
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%             (optional) 
    

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/22/2017        0.9             hme
%                


v=clv.tuvalues;
N=clv.tusize;
vi=clv.tuvi;

if nargin<2
    error('A game and a coalition structure cs must be provided!')
elseif nargin < 3
  x=clv.weightedKernel(cs);
  warning('PKQ:NoPayoffInput','Computing default payoff!');
  n=length(x);
  tol=10^6*eps;
  pS='';
elseif nargin< 4
   n=length(x);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
  tol=10^6*eps;
  pS='';
elseif nargin< 5
   n=length(x);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
  tol=10^6*eps;
else
   n=length(x);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
end
smat=-inf;
if iscell(cs)
   cs=clToMatlab(cs);
else
  cs=double(cs);
end
lcs=length(cs);

if isempty(pS)
   S=1:N;
   for k=1:n, A1(:,k) = -bitget(S,k);end
   warning('No vector of weights specified assuming per capita pre-kernel.')
   mat=-A1';
   clS=ones(1,n)*mat;
   pS=1./clS;
   pS(N)=1;
end

[smat,ei,ecs]=cs_msrpls(v,x,n,pS,cs);
effQ=all(abs(ecs)<tol);
krQ=false;
if effQ==0, return; end

%
% Selecting only the l,k in B from smat.

mb=1:n;
lcs=length(cs);
for kk=1:lcs
    clm{kk}=mb(bitget(cs(kk),mb)==1);
end

cssm=diag(ei);

for kk=1:lcs
    ckk=clm{kk};
    lckk=numel(clm{kk});
    if lckk>1
       for ii=1:lckk-1
           for jj=ii+1:lckk
               cssm(ckk(ii),ckk(jj))=smat(ckk(ii),ckk(jj));
               cssm(ckk(jj),ckk(ii))=smat(ckk(jj),ckk(ii));
           end
       end
    end
end
smat=cssm;


ir=x-vi;
irQ=all(ir>=-tol);
if irQ
   krm=smat-smat';
   irm=repmat(ir,n,1);
   kriQ=all((krm.*irm)<=tol);
   krQ=all(kriQ);
else
  kriQ=0;
  krQ=0;
end


%-------------------------------------
function [smat,ei,ecs]=cs_msrpls(v,x,n,pS,cs) 
% Computes the maximum surpluses w.r.t. payoff x.
% output:
%  smat     -- Matrix of maximum surpluses.
%  ei       -- The excess of singleton coalitions.       
%  ecs      -- The excess of the coaltions belonging to cs.
%
%  input:
%  v      -- A Tu-Game v of length 2^n-1.
%  x      -- payoff vector of length(1,n)
%  n      -- number of players
%  pS     -- A vector of weights of length 2^n-1.

% the excesses of x wrt. the game v
% Borrowed from J. Derks
Xm=x(1); for ii=2:n, Xm=[Xm x(ii) Xm+x(ii)]; end
e=v-Xm;
pl=1:n;
ci=2.^(pl-1);
e=v-Xm;
ei=e(ci);
ecs=e(cs);
clear v Xm;
e=e.*pS;
e(cs)=ecs;
% Determing max surpluses
[se, sC]=sort(e,'descend');
smat=-inf(n);
q0=n^2-n;
q=0;
k=1;
while q~=q0
  kS=sC(k);
  ai=bitget(kS,pl)==1;
  bj=ai==0;
  pli=pl(ai);
  plj=pl(bj);
  if isempty(plj)==0
    for i=1:numel(pli)
      for j=1:numel(plj)
        if smat(pli(i),plj(j))==-Inf
           smat(pli(i),plj(j))=se(k); % max surplus of i against j.
           q=q+1;
        end
      end
    end
  end
  k=k+1;
end
smat=tril(smat,-1)+triu(smat,1);

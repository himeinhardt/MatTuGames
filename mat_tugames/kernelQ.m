function [krQ, kriQ, smat]=kernelQ(v,x,tol)
% KERNELQ checks whether the imputation x is a kernel element 
% of the TU-game v.
% 
%  Usage:[krQ, kriQ, smat]=kernelQ(v,x,tol);
%
%
% Define variables:
%  output:
%  krQ      -- Returns 1 (true) whenever the impuatation x is 
%              a kernel element, otherwise 0 (false).
%  kriQ     -- Returns a list of true (1) and/or false (0) of length n.
%  smat     -- Matrix of maximum surpluses.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%             (optional) 
    

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/17/2013        0.3             hme
%                

if nargin < 2
  x=PreKernel(v);
  warning('PKQ:NoPayoffInput','Computing default payoff!');
  n=length(x);
  tol=10^6*eps;
elseif nargin< 3
   N=length(v);
   n=length(x);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
  tol=10^6*eps;
else
   N=length(v);
   n=length(x);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
end
smat=-inf;
e=0;
effQ=abs(v(end)-sum(x))<tol;
krQ=0;
if effQ==0, return; end 

smat=msrpls(v,x,n);
l=1:n;
ir=x-v(bitset(0,l));
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
function smat=msrpls(v,x,n) 
% Computes the maximum surpluses w.r.t. payoff x.
% output:
%  smat     -- Matrix of maximum surpluses.
%
%  input:
%  v      -- A Tu-Game v of length 2^n-1.
%  x      -- payoff vector of length(1,n)


% the excesses of x wrt. the game v
% Borrowed from J. Derks
Xm=x(1); for ii=2:n, Xm=[Xm x(ii) Xm+x(ii)]; end
e=v-Xm;
clear v Xm;
% Determing max surpluses
[se, sC]=sort(e,'descend');
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

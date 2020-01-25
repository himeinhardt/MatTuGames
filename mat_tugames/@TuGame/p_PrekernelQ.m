function [prkQ, smat, e]=p_PrekernelQ(clv,x,tol)
% P_PREKERNELQ checks whether the imputation x is a pre-kernel element 
% of the TU-game v.
%
%  Usage: [prkQ smat e]=p_PrekernelQ(clv,x,tol); 
%
%
% Define variables:
%  output:
%  prkQ     -- Returns 1 (true) whenever the impuatation x is 
%              a pre-kernel element, otherwise 0 (false).
%  smat     -- Matrix of maximal surpluses.
%  input:
%  clv      -- TuGame class object.
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
%   10/30/2012        0.3              hme
%                

if nargin < 2 
   if isa(clv,'TuSol')
      x=clv.tu_prk;
   elseif isa(clv,'p_TuSol')
      x=clv.tu_prk;
   else
      x=clv.p_PreKernel();
   end
   if isempty(x)
     x=clv.p_PreKernel();
   end
   warning('PPKQ:NoPayoffInput','Computing/Using default payoff!');
   tol=10^6*eps;
elseif nargin<3
  tol=10^6*eps;
else
  tol=10^6*eps;
end

v=clv.tuvalues;
n=clv.tuplayers;

smat=-inf;
e=0;
effQ=abs(v(end)-sum(x))<tol;
prkQ=false;
if effQ==0, return; end

[e, smat]=p_msrpls(v,x,n);
e=e';
lms=abs(smat-smat')<tol;
prkQ=all(all(lms));


%-------------------------------------
function [e, smat]=p_msrpls(v,x,n)
% Computes the maximum surpluses w.r.t. payoff x.
% output:
%  e        -- Excess vector
%  smat     -- Matrix of maximal surpluses.
%
%  input:
%  v      -- A Tu-Game v of length 2^n-1.
%  x      -- payoff vector of length(1,n)

% Borrowed from J. Derks
Xm{1}=x(1); for ii=2:n, Xm{1}=[Xm{1} x(ii) Xm{1}+x(ii)]; end
% Determining max surpluses.
spmd
 cv=codistributed(v);
 cXm=codistributed(Xm{1});
 e=cv-cXm;
 [se1, sC1]=sort(e,'descend');
 Smat=-inf(n);
end
sC=gather(sC1);
se=gather(se1);
q0=n^2-n;
q=0;
k=1;
pl=1:n;
spmd
 while q~=q0
   kS=sC(k);
   ai=bitget(kS,pl)==1;
   bj=ai==0;
   pli=pl(ai);
   plj=pl(bj);
   if isempty(plj)==0
     for i=1:numel(pli)
       for j=1:numel(plj)
         if Smat(pli(i),plj(j))==-Inf
            Smat(pli(i),plj(j))=se(k); % max surplus of i against j.
            q=q+1;
         end
       end
     end
   end
   k=k+1;
 end
end
smat=Smat{1};
smat=tril(smat,-1)+triu(smat,1);

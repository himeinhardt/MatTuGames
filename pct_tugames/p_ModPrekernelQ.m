function [prkQ smat e]=p_ModPrekernelQ(v,x,tol)
% P_MODPREKERNELQ checks whether the imputation x is a modified pre-kernel element 
% of the TU-game v using Matlab's PCT.
% 
%  Usage:[prkQ smat e]=p_ModPrekernelQ(v,x,tol);
%
%
% Define variables:
%  output:
%  prkQ     -- Returns 1 (true) whenever the impuatation x is 
%              a modified pre-kernel element, otherwise 0 (false).
%  smat     -- Matrix of maximum surpluses.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%             (optional) 
%
    

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/11/2018        1.0              hme
%                


if nargin < 2
  x=p_PreKernel(v);
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
prkQ=false;
if effQ==0, return; end


[e, smat]=p_msrpls(v,x);
e=e';
lms=abs(smat-smat')<tol;
prkQ=all(all(lms));


%-------------------------------------
function [e, smat]=p_msrpls(v,x)
% Computes the maximum surpluses w.r.t. payoff x.
% output:
%  e        -- Excess vector
%  smat     -- Matrix of maximum surpluses.
%
%  input:
%  v      -- A Tu-Game v of length 2^n-1.
%  x      -- payoff vector of length(1,n)


N=length(v); n=length(x);
dv=dual_game(v);
% the excesses of x wrt. the game v
% Borrowed from J. Derks
Xm{1}=x(1); for ii=2:n, Xm{1}=[Xm{1} x(ii) Xm{1}+x(ii)]; end
% Determining max surpluses.
v=parallel.pool.Constant(v);
dv=parallel.pool.Constant(dv);
cXm=parallel.pool.Constant(Xm{1});
clear Xm;
spmd
%cv=codistributed(v);
%cXm=codistributed(Xm{1});
e=v.Value-cXm.Value;
de=dv.Value-cXm.Value;
eh = max(e);
deh = max(de);
e=max(e+deh,de+eh);
[se1, sC1]=sort(e,'descend');
Smat=-inf(n);
end
sC=gather(sC1);
se=gather(se1);
clear sC1 se1;
sC=parallel.pool.Constant(sC);
se=parallel.pool.Constant(se);
q0=n^2-n;
q=0;
k=1;
pl=1:n;
spmd
 while q~=q0
   kS=sC.Value(k);
   ai=bitget(kS,pl)==1;
   bj=ai==0;
   pli=pl(ai);
   plj=pl(bj);
   if isempty(plj)==0
     for i=1:numel(pli)
       for j=1:numel(plj)
         if Smat(pli(i),plj(j))==-Inf
            Smat(pli(i),plj(j))=se.Value(k); % max surplus of i against j.
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

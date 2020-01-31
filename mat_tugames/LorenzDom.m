function LD = LorenzDom(v,x,y,tol)
% LORENZDOM checks if x Lorenz dominates y in game v.
% 
%
% Usage: [crq x]=LorenzDom(v,x,y,tol)
%
% Define variables:
%  output:
%  crq      -- Returns 1 (true) or 0 (false).
%  x        -- The smallest allocation that satisfies all core constraints.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/12/2015        0.7             hme
%   01/15/2019        1.0             hme
%



if nargin < 2
   error('Game must be given!');
elseif nargin< 3
   N=length(v);
   n=length(x);
   msg='Using the Shapley value as second payoff.';
   warning('LorD:SecPayoff',msg);
   y=ShapleyValue(v);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
  tol=10^6*eps;
elseif nargin< 4
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

ldQ=false;
dij.x=DIJ(v,x,n);
dij.y=DIJ(v,y,n);


phi.x=phif(sort(x),n,N);
phi.y=phif(sort(y),n,N);
ldQ=all(phi.x>=phi.y-tol);
LD=struct('ldQ',ldQ,'phi',phi,'dij',dij);


%----------------------
function dij=DIJ(v,z,n)

[~,smat]=PrekernelQ(v,z);
%dij=zeros(n,n);
del=zeros(n,n);

for ii=1:n-1
   for jj=ii+1:n
           del(ii,jj)=min((z(jj)-z(ii))/2,-smat(ii,jj));
           dij(ii,jj)=max(del(ii,jj),0);
           del(jj,ii)=min((z(ii)-z(jj))/2,-smat(jj,ii));
           dij(jj,ii)=max(del(jj,ii),0);            
   end
end



%----------------------
function phi=phif(x,n,N)
S=1:N;
it=0:-1:1-n;
sS=rem(floor(S(:)*pow2(it)),2);
ov=ones(n,1);
csz=sS*ov;
phi=zeros(1,n);
adx=additive_game(x);
for k=1:n 
    kS=S(csz==k);
    phi(k)=min(adx(kS));
end

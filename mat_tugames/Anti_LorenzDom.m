function ALD = Anti_LorenzDom(v,x,y,tol)
% Anti_LORENZDOM checks if x anti-Lorenz dominates y in game v.
%
% Source: Hougaard, J.L., Peleg, B., Thorlund-Petersen, L., 2001. On the set of Lorenz-maximal
%         imputations in the core of a balanced game. Internat. J. Game Theory 30, 147–165.
%
% Usage: ALD=Anti_LorenzDom(v,x,y,tol)
%
% Define variables:
%  output: field variables
%  ldQ      -- Returns 1 (true) if x Lorenz Dom y, otherwise 0 (false).
%              Not required that y is in the core.
%  phi      -- Mapping of phi(x) and phi(y). 
%  dij      -- Matrices of dij(x)=max(del_ij(x),0) and dij(y)=max(del_ij(y),0) 
%              for all pairs (i,j). Zero matrix is a necessary condition for x to 
%              be Lorenz max within the core, whenever x is in the core.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- An allocation x of length n.
%  y        -- An alternative allocation of length n.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/20/2021        1.9.1           hme
%



if nargin < 2
   error('Game must be given!');
elseif nargin< 3
   N=length(v);
   n=length(x);
   msg='Using the Shapley value as second payoff.';
   warning('ALorD:SecPayoff',msg);
   y=ShapleyValue(v);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
  tol=10^8*eps;
elseif nargin< 4
   N=length(v);
   n=length(x);
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
  tol=10^8*eps;
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

phi.x=phif(x,n,N);
phi.y=phif(y,n,N);
ldQ=all(phi.x<=phi.y+tol);
ALD=struct('ldQ',ldQ,'phi',phi,'dij',dij);


%----------------------
function dij=DIJ(v,z,n)

[~,smat]=Anti_PrekernelQ(v,z);
dij=zeros(n,n);
del=zeros(n,n);

if n>1
  for ii=1:n-1
      for jj=ii+1:n
           del(ii,jj)=max((z(ii)-z(jj))/2,smat(jj,ii));
           dij(ii,jj)=min(del(ii,jj),0);
           del(jj,ii)=max((z(jj)-z(ii))/2,smat(ii,jj));
           dij(jj,ii)=min(del(jj,ii),0);            
      end
  end
else
  dij=min(smat,0)
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
    phi(k)=max(adx(kS));
end

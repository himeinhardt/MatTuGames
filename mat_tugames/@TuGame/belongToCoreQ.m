function [crq, bS]=belongToCoreQ(clv,x,str,tol)
% BELONGTOCOREQ checks whether the imputation x belongs to the core of game v.
% 
%
%  Usage: [crq bS]=clv.belongToCoreQ(x,str,tol)
%
% Define variables:
%  output:
%  crq      -- It returns 1 (true) or 0 (false).
%  bS       -- List of blocking coalitions.
%
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n) 
%  str      -- This string defines rational approximation.
%              Permissible methods are:
%              'rat' for a rational approximation
%              'real' otherwise.
%              Default is 'rat'.
%  tol      -- A positive tolerance value. The default is set to
%              (2*10^4*eps) for method 'rat',
%              (10^8*eps) for mehtod 'real'.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

S=1:N;
if nargin<2
  if isa(clv,'TuSol')
     x=clv.tu_prk;
  elseif isa(clv,'TuVal')
     x=clv.tu_sh;
  elseif isa(clv,'p_TuVal')
     x=clv.tu_sh;
  elseif isa(clv,'TuGame')
     x=clv.PreKernel();
  end
  impv=num2str(x);
  msg='A second input argument was not given. Using default imputation: \n\n %s \n';
  warning(msg,impv) ;
  str='rat'; % assuming rational numbers.
  tol=2*10^4*eps;
elseif nargin<3
   str='rat'; % assuming rational numbers.
   tol=2*10^4*eps;
elseif nargin<4
   tol=10^8*eps;
else
end
if ischar(x) % Convert a rational number into reals.
    x=str2num(x);
end
tol=-tol;
grq=false(N,1);
crq=false(1,N);
sumx=x(1); for ii=2:n, sumx=[sumx x(ii) sumx+x(ii)]; end

if strcmp(str,'rat')
 grq=sumx>=v+tol;
 grq(N)=abs(sumx(N)-v(N))<=abs(tol);
elseif strcmp(str,'real')
 dsv=sumx-v;
 grq=dsv>tol;
 grq(N)=abs(sumx(N)-v(N))<=abs(tol);
else
 dsv=sumx-v;
 grq=dsv>tol;
 grq(N)=abs(sumx(N)-v(N))<=abs(tol);
end

crq=all(grq');
lg=grq==0;
bS=S(lg);

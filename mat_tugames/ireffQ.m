function IREFF=ireffQ(v,x,tol)
% IREFFQ checks if a payoff vector satisfies the individual rationalitiy as well as
% the efficiency property.
%
%  USAGE: IREFF=ireffQ(v,x)
%
%
% Define variables:
%  output: Fields
%  irQ      -- Returns true (1), if x is individual rational, otherwise
%              false (0).
%  effQ     -- Returns true (1), if x is efficient, otherwise
%              false (0).
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 
%              




%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/22/2017        0.9             hme
%



if nargin<3
 tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
k=1:n;
ic=2.^(k-1);
vi=v(ic);
irQ=all(x+tol>=vi);
effQ=abs(sum(x)-v(N))<tol;

IREFF=struct('irQ',irQ,'effQ',effQ);


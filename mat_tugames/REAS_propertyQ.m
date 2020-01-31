function REAS=REAS_propertyQ(v,x,tol)
% REAS_PROPERYQ checks if the vector x satisfies the reasonableness
% on both sides (REAS).
%
% Usage: REAS=REAS_propertyQ(v,x,tol)
%
% Define structure variables:
%  output:
%  reasQ    -- Returns true (1) if the solution x satisfies REAS,
%              otherwise false (0).  
%  ub       -- REAS from above.
%  lb       -- REAS from below. 
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/10/2017        1.0             hme
% 


if nargin < 3
    tol=10^8*eps;
end    

mgc=AllMarginalContributions(v);
c=max(mgc);
f=min(mgc);

rabQ=all(x<=c+tol);
rfbQ=(all(f<=x+tol));

REAS.reasQ=(rabQ && rfbQ);
REAS.ub=c;
REAS.lb=f;

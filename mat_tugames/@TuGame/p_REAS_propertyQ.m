function REAS=p_REAS_propertyQ(clv,x,tol)
% P_REAS_PROPERYQ checks if the vector x satisfies the reasonableness
% on both sides (REAS) using Matlab's PCT.
%
% Usage: REAS=clv.p_REAS_propertyQ(x,tol)
%
% Define structure variables:
%  output:
%  reasQ      -- Returns true (1) if the solution x satisfies REAS,
%                otherwise false (0).   
%  input:
%  clv      -- TuGame class object.
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

mgc=clv.p_AllMarginalContributions();
c=max(mgc);
f=min(mgc);

rabQ=all(x<=c+tol);
rfbQ=(all(f<=x+tol));

REAS.reasQ=(rabQ && rfbQ);


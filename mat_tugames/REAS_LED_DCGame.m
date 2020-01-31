function [LED,REAS]=REAS_LED_DCGame(v,x,tol)
% REAS_LED_DCGAME verifies that a reasonable vector x of game v, then the shifted 
% ducal cover game statisfies LED w.r.t. the replicated vector (x,x).
%
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          Sudhoelter (1997), The modified nucleolus: Properties and axiomatizations. International Journal of Game Theory, 26
%          (2):147–182, Jun 1997. ISSN 1432-1270. doi: 10.1007/BF01295846. URL https://doi.org/10.1007/BF01295846.

%
% Usage: [LED,REAS]=REAS_LED_DCGame(v,x,tol)
%
% Define variables:
%
%  Output structure variables:
% LED
%  propQ    -- Returns true (1) whenever the solution has large excess difference,
%              otherwise false (0).
%  lexd     -- Retruns the value of the large excess difference.
%  x        -- Replicates the input variable x.
%
% REAS
%  reasQ     -- Returns true (1) if the solution x satisfies REAS,
%                otherwise false (0).  
%  ub       -- REAS from above.
%  lb       -- REAS from below
%
%  Input:
%  v        -- A Tu-Game v of length 2^n-1.
%  x        -- payoff vector of size(1,n) (optional)
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/15/2018        1.0             hme



if nargin<3
 tol=10^6*eps;
elseif nargin==3
 if ischar(tol)
    error('Tolerance value is a string character! Must be a number!');
 end
end

REAS=REAS_propertyQ(v,x,tol);


if REAS.reasQ==1
   dc_v=DualCover(v);
   if length(v)>1
      Mgc=AllMarginalContributions(v);
      d=max(Mgc)-min(Mgc);
      t=6*sum(d);
   else
      d=0;
      t=0;
   end
end

v_x=shiftGame(dc_v,t);
y=[x,x];
LED=LED_propertyQ(v_x,y,tol);

function tB=p_totallyBalancedQ(clv,tol)
% P_TOTALLYBALANCEDQ checks whether the core of all subgames is non-empty
% using Matlab's PCT.
%
%
% Usage: tB=clv.p_totallyBalancedQ(tol)
% Define structure variables:
%  output:
%  Q         -- The game is totally balanced if Q=1, otherwise false.
%  sbgQ       -- An array of ones (true) and/or false indicating if the 
%               the core of the associated subgame is empty or not. 
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/28/2018        0.9             hme
%

if nargin<2
 tol=10^10*eps; % Change this value if the solution is not correct.
end

N=clv.tusize;
crQ=false(1,N);

parfor S=1:N
    subg=clv.SubGame(S);
    try  
      crQ(S)=CddCoreQ(subg);
    catch
      crQ(S)=coreQ(subg,tol);
    end
end

tB.Q=all(crQ);
tB.sbgQ=crQ;

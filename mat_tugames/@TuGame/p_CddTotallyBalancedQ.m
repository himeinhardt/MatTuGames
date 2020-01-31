function tB=p_CddTotallyBalancedQ(clv,tol)
% P_CDDTOTALLYBALANCEDQ checks whether the core of all subgames is non-empty using cddmex
% and Matlab's PCT.
%
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% and the Matlab interface
% to the cdd solver (cddmex) http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Usage: tB=clv.p_CddTotallyBalancedQ(tol)
% Define structure variables:
%  output:
%  Q         -- The game is totally balanced if Q=1, otherwise false.
%  crQ       -- An array of ones (true) and/or false indicating if the 
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
 tol=10^6*eps; % Change this value if the solution is not correct.
end

N=clv.tusize;
crQ=false(1,N);

parfor S=1:N
    subg=clv.SubGame(S);
    crQ(S)=CddCoreQ(subg,tol);
end

tB.Q=all(crQ);
tB.sbgQ=crQ;

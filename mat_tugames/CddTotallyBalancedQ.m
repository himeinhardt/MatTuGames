function tB=CddTotallyBalancedQ(v,tol)
% CDDTOTALLYBALANCEDQ checks whether the core of all subgames is non-empty using cddmex.
%
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% and the Matlab interface
% to the cdd solver (cddmex) http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Usage: tB=CddTotallyBalancedQ(v,tol)
% Define structure variables:
%  output:
%  Q         -- The game is totally balanced if Q=1, otherwise false.
%  crQ       -- An array of ones (true) and/or false indicating if the 
%               the core of the associated subgame is empty or not. 
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
%   01/28/2018        0.9             hme
%

if nargin<2
 tol=10^6*eps; % Change this value if the solution is not correct.
end

N=length(v);
%[~, n]=log2(N);
crQ=false(1,N);

for S=1:N
    subg=SubGame(v,S);
    crQ(S)=CddCoreQ(subg,tol);
end

tB.Q=all(crQ);
tB.sbgQ=crQ;

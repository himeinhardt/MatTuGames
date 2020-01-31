function tB=p_totallyBalancedQ(v,tol)
% P_TOTALLYBALANCEDQ checks whether the core of all subgames is non-empty
% using Matlab's PCT.
%
%
% Usage: tB=p_totallyBalancedQ(v,tol)
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
 tol=10^10*eps; % Change this value if the solution is not correct.
end

N=length(v);
%[~, n]=log2(N);
crQ=false(1,N);

parfor S=1:N
    subg=SubGame(v,S);
    try  
      crQ(S)=CddCoreQ(subg);
    catch
      crQ(S)=coreQ(subg,tol);
    end
end

tB.Q=all(crQ);
tB.sbgQ=crQ;

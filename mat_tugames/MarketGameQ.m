function mgQ=MarketGameQ(v,tol)
% MARKETGAMEQ checks whether the game v is a market game.
%
%
% Usage: mgQ=MarketGameQ(v,tol)
% Define structure variables:
%  output:
%  Q         -- Returns true (1), if the game is a market game, i.e., is totally balanced, otherwise false.
%  crQ       -- An array of ones (true) and/or false indicating if the 
%               the core of the associated subgame is empty or not.
%
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/24/2022        1.9.1           hme
%    
msg=nargchk(1,2,nargin);
error(msg);

if nargin<2
 tol=10^8*eps;
end    
mgQ=totallyBalancedQ(v,tol);    

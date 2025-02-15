function DPG=DecomposeInPositiveGames(v,tol)
% DecomposeInPositiveGames decomposes a TU game v
% into the difference of two positive games (convex games).
%
% Usage: DPG=DecomposeInPositiveGames(v)
%
% Define variables:
%  output:  structure elements
%  eqQ      -- Returns true (1), if the original and the recomposed 
%              game are equal.
%  dp1      -- First decomposed positive game.    
%  dp2      -- Second decomposed positive game.    
%  w        -- Recomposed game.
%  v        -- Original game.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
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
%   02/24/2022        1.9.1           hme
%
     
if nargin < 2
   tol=10^6*eps; 
end    
    
[hd utmat]=unanimity_games(v);
lhd=max(0,hd);
shd=-min(0,hd);
pg1=utmat*lhd';
pg2=utmat*shd';
w=(pg1-pg2)';
eqQ=all(abs(v-w)<tol);
DPG.eqQ=eqQ;
DPG.dp1=pg1';
DPG.dp2=pg2';
DPG.v=v;
DPG.w=w;


end

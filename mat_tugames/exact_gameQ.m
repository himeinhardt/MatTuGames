function egQ=exact_gameQ(v,tol);
% EXACT_GAMEQ checks whether game v is an exact game using 
% the Matlab's Optimization toolbox. Uses Dual-Simplex (Matlab R2015a).
% 
%  Usage: [v_e xm status]=exact_gameQ(v,tol);
%
% Define variables:
%  output:
%  egQ      -- Returns one if the game is exact, otherwise zeros.
%    
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to -10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/01/2018        1.0             hme
%

if nargin<2
 tol=10^6*eps;
end

v_e=exact_game(v,tol);
egQ=all(abs(v-v_e)<tol);

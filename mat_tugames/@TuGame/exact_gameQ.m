function egQ=exact_gameQ(clv,tol);
% EXACT_GAMEQ checks whether game v is an exact game using 
% the Matlab's Optimization toolbox. Uses Dual-Simplex (Matlab R2015a).
% 
%  Usage: egQ=clv.exact_gameQ(tol);
%
% Define variables:
%  output:
%  egQ      -- Returns one if the game is exact, otherwise zeros.
%    
%  input:
%  clv      -- TuGame class object.
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
v=clv.tuvalues;
v_e=clv.exact_game(tol);
egQ=all(abs(v-v_e)<tol);

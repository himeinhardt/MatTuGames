function hd=p_harsanyi_dividends(v)
% P_Harsanyi_dividends computes the unanimity coordinates or Harsanyi dividends.
% Same as UNANIMITY_GAMES. Call the function unanimity_games directly.
% Using Matlab's PCT.
%
% Usage: hd=p_harsanyi_dividends(v)
% Define variables:
%  output:
%  hd       -- Unanimity coordinates or Harsanyi dividends.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/21/2011        0.1 alpha        hme
%                


hd=p_unanimity_games(v);

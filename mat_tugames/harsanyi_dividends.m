function hd=harsanyi_dividends(v)
% Harsanyi_dividends computes the unanimity coordinates or Harsanyi dividends.
% Same as UNANIMITY_GAMES. Call the function unanimity_games directly.
% For n>14 this function needs some time to complete.
%
% Usage: hd=harsanyi_dividends(v)
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
%   08/09/2010        0.1 beta        hme
%                


hd=unanimity_games(v);

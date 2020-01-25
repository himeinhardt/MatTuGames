function hd=harsanyi_dividends(clv)
% Harsanyi_dividends computes the unanimity coordinates or Harsanyi dividends.
% Same as UNANIMITY_GAMES. Call the function unanimity_games directly.
% For n>14 this function needs some time to complete.
%
% Usage: hd=harsanyi_dividends(clv)
%
% Define variables:
%  output:
%  hd       -- Unanimity coordinates or Harsanyi dividends.
%
%  input:
%  clv        -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%                


hd=unanimity_games(clv);

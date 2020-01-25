function sv=simple_game(w_coal,n)
% SIMPLE_GAME computes from a list of winning coalitions the 
% corresponding simple game.
%
% Usage: sv=simple_game(w_coal,n)
% Example:
% Let w_coal=[10 12 7 11 13 14 15] be a set of winning coalitions, then
% sv=simple_game(w_coal,4);
% computes the corresponding simple game.
%
%
% Define variables:
%  output:
%  sv       -- A simple game.
%  input:
%  w_coal   -- A list/vector of winning coalitions
%  n        -- The number of players involved (positive number).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/13/2010        0.1 beta        hme
%   05/15/2014        0.5             hme
%                


S=1:2^n-1;
w_coal=sort(w_coal);
sv=ismembc(S,w_coal);

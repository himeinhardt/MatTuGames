function [v,hd]=product_game(x)
% PRODUCT_GAME computes from a vector x the corresponding product game.
%
% Source: Rosales, D. (2014), Cooperative product games (www.academia.edu/401764) 
%         Pierre Dehez (2023), COOPERATIVE PRODUCT GAMES, LIDAM Discussion Paper CORE 2023 / 10 
%
% Usage: v=product_game(x)
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1.
%  hd       -- Unanimity coordinates or Harsanyi dividends.
%
% input:
%  x        -- a vector of size(1,n). (players weights)
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/05/2015        0.7               hme
%   04/04/2021        1.9.1             hme
%

n=length(x);
v=x(1); for ii=2:n, v=[v x(ii) v*x(ii)]; end
v=v-1;     % Proposed adjustment by P. Dehez.  
hd=x(1)-1; for ii=2:n, hd=[hd x(ii)-1 hd*(x(ii)-1)]; end

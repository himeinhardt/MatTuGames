function v=product_game(x)
% PRODUCT_GAME computes form a vector x the corresponding product game.
%
% Usage: v=product_game(x)
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1.
%
% input:
%  x        -- a vector of size(1,n).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/05/2015        0.7             hme
%

n=length(x);
v=x(1); for ii=2:n, v=[v x(ii) v*x(ii)]; end

function v=additive_game(x)
% ADDITIVE_GAME computes form a vector x the corresponding additive game.
%
% Usage: v=additive_game(ad)
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
%   11/23/2012        0.2             hme
%

n=length(x);
% Computing the additive game w.r.t. x.
% Borrowed from J. Derks
v=x(1); for ii=2:n, v=[v x(ii) v+x(ii)]; end

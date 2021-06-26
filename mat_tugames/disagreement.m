function da=disagreement(v)
% DISAGREEMENT computes the disagreement vector of game v.
%
% Usage: da=disagreement(v)
% Define variables:
%  output:
%  da       -- The disagreement vector of game v.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08.27.2020        1.9             hme
%

[~,da]=UtopiaPayoff(v);

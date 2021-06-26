function da=disagreement(clv)
% DISAGREEMENT computes the minimum claim vector of game v.
%
% Usage: da=clv.disagreement()
% Define variables:
%  output:
%  da       -- The disagreement vector of game v.
%
%  input:
%  clv      -- TuGame class object.
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

[~,da]=clv.UtopiaPayoff();

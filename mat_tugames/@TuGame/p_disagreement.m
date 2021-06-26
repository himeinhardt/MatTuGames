function da=p_disagreement(clv)
% P_DISAGREEMENT computes the disagreement vector of game v using MATLAB's PCT.
%
% Usage: da=clv.p_disagreement()
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

[~,da]=clv.p_UtopiaPayoff();

function v=probability_game(x)
%PROBABILITY_GAME computes from a vector x the corresponding probability game.
%
% Source: Hou, D., Xu, G., Sun, P., Driessen, T., 2018. The Shapley value for the probability game. Oper. Res. Lett. 46, 457–461.
%         Pierre Dehez (2023), Sharing a collective probability of success, Mathematical Social Sciences 123, 122–127.
%
% Usage: v=probability_game(x)
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1.
%
% input:
%  x        -- a vector of size(1,n) between zero and one. (players weights)
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/20/2021        1.9.1             hme
%

n=length(x);
v=1-x(1); for ii=2:n, v=[v 1-x(ii) v*(1-x(ii))]; end
v=1-v;  

function [vR,mnbp]=root_game(clv)
% ROOT_GAME computes from game v its associated root game vR.
%
% Usage: clv=root_game()
%
%  Source: Zhao (2001) and Calleja et. al. (2009).
%
% Define variables:
%  output:
%  vR        -- A root game of length 2^n-1.
%  mnbp      -- The minimum no blocking payoff.
%
%  input:
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/24/2017        0.9             hme
%                

v=clv.tuvalues;
N=clv.tusize;
MNBP_vR=clv.minNoBlockPayoff();
vR=v;
mnbp=MNBP_vR.mnbp;
vR(N)=mnbp;

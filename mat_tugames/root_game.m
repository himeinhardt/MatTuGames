function [vR,mnbp]=root_game(v)
% ROOT_GAME computes from game v its associated root game vR.
%
% Usage: v=root_game(v)
%
%  Source: Zhao (2001) and Calleja et. al. (2009).
%
% Define variables:
%  output:
%  vR        -- A root game of length 2^n-1.
%  mnbp      -- The minimum no blocking payoff.
%  input:
%  v         -- A TU game v of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/22/2017        0.9             hme
%                

N=length(v);
MNBP_vR=minNoBlockPayoff(v);
vR=v;
mnbp=MNBP_vR.mnbp;
vR(N)=mnbp;

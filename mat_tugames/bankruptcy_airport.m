function [bv,av]=bankruptcy_airport(E,d_vec)
% BANKRUPTCY_AIRPORT computes from a bankruptcy problem the associated 
% bankruptcy as well as the airport profit game. Both games must
% be equal.
%
% Usage: [bv,av]=airport_profit(E,d_vec)
%
%
% Define variables:
%  output:
%  bv        -- A TU-bankruptcy game.
%  av        -- The airport surplus game v of length 2^n-1.
%
%  input:
%  c        -- A cost vector of length n.
%  b        -- A benefit vector of length n.
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/26/2015        0.6             hme
%

n=length(d_vec);
Ci=sum(d_vec)-E;
Ci=ones(1,n)*Ci;
bi=d_vec;
av=airport_profit(Ci,bi);
bv=bankruptcy_game(E,d_vec);

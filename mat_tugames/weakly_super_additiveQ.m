function wsaq=weakly_super_additiveQ(v)
% WEAKLY_SUPER_ADDITIVEQ returns 1 whenever the game v is weakly super
% additive. 
%
% wsaq=weakly_super_additiveQ(v)
% Define variables:
%  output:
%  wsaq     -- Returns 1 (true) or 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/06/2010        0.1 beta        hme
%                


dv=dual_game(v);
lv=v<=dv;
wsaq=all(lv);

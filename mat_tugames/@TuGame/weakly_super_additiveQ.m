function wsaq=weakly_super_additiveQ(clv)
% WEAKLY_SUPER_ADDITIVEQ returns 1 whenever the game v is weakly super
% additive. 
%
% wsaq=weakly_super_additiveQ(clv)
%
% Define variables:
%  output:
%  wsaq     -- Returns 1 (true) or 0 (false).
%
%  input:
%  clv        -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%                

v=clv.tuvalues;
dv=dual_game(clv);
lv=v<=dv;
wsaq=all(lv);

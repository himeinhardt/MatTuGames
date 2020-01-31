function wsaq=weakly_sub_additiveQ(v,tol)
% WEAKLY_SUB_ADDITIVEQ returns 1 whenever the game v is weakly 
% sub-additive. 
%
% wsaq=weakly_sub_additiveQ(v)
% Define variables:
%  output:
%  wsaq     -- Returns 1 (true) or 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- Tolerance value. By default, it is set to 2*10^4*eps.
%              (optional)
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/20/2019        1.1             hme
%                

if nargin<2
   tol=2*10^4*eps;
end


dv=dual_game(v);
lv=v+tol>=dv;
wsaq=all(lv);

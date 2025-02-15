function [rtQ,v_r]=rootedQ(v,tol)
% ROOTEDQ checks whether the game v is rooted.
%
%
% Usage: rtQ=rootedQ(v)
% Define variables:
%  output:
%  rt       -- Returns 1 (true) whenever the game v is rooted,
%                otherwise 0 (false).
%  v_r      -- The unique rooted TU-game of v.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- Tolerance value. By default, it is set to 10^4*eps.
%              (optional)
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/06/2022        1.9.1           hme
%


if nargin < 2
  tol=10^4*eps;
end	

v_r=root_game(v);
rtQ=all(abs(v-v_r)<=tol);

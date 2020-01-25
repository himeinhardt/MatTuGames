function zmq=p_zero_monotonicQ(v);
% P_ZERO_MONOTONICQ returns 1 whenever the game v is zero-monotonic.
% Using Matlab's PCT.
%
% Usage: zmq=p_zero_monotonicQ(v)
% Define variables:
%  output:
%  zmq      -- Returns 1 (true) or 0 (false).
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
%   05/20/2011        0.1 alpha        hme
%                


zv=zero_normalization(v);
zmq=p_monotone_gameQ(zv);

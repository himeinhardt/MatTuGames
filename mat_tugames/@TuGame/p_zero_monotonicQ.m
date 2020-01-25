function zmq=p_zero_monotonicQ(clv);
% P_ZERO_MONOTONICQ returns 1 whenever the game v is zero-monotonic.
% Using Matlab's PCT.
%
% Usage: zmq=p_zero_monotonicQ(clv)
%
% Define variables:
%  output:
%  zmq      -- Returns 1 (true) or 0 (false).
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
%   10/28/2012        0.3             hme
%                

v=clv.tuvalues;

zv=zero_normalization(v);
zmq=p_monotone_gameQ(zv);

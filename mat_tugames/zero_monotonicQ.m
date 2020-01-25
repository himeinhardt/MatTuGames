function zmq=zero_monotonicQ(v,tol);
% ZERO_MONOTONICQ returns 1 whenever the game v is zero-monotonic.
%
% Usage: zmq=zero_monotonicQ(v)
% Define variables:
%  output:
%  zmq      -- Returns 1 (true) or 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- Tolerance value. By default, it is set to -2*10^4*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/13/2010        0.1 beta        hme
%   06/01/2016        0.8             hme
%                

if nargin<2
   tol=-2*10^4*eps;
end


zv=zero_normalization(v);
zmq=monotone_gameQ(zv,tol);

function zmq=zero_monotonicQ(clv,tol);
% ZERO_MONOTONICQ returns 1 whenever the game v is zero-monotonic.
%
% Usage: zmq=zero_monotonicQ(clv)
%
% Define variables:
%  output:
%  zmq      -- Returns 1 (true) or 0 (false).
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. By default, it is set to -2*10^4*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/28/2012        0.3             hme
%   06/01/2016        0.8             hme
%                

if nargin<2
   tol=-2*10^4*eps;
end


zv=clv.zero_normalization;
zmq=monotone_gameQ(zv);

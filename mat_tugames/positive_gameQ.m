function psQ=positive_gameQ(v,tol);
% POSITIVE_GAMEQ returns true (1) if all Harsanyi dividends are non-negative, otherwise false (zero). 
%
% Usage: psQ=positive_gameQ(v,tol);
% Define variables:
%  output:
%  psQ      -- Returns 1 (true) or 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- Tolerance value. By default, it is set to (10^6*eps).
%              (optional) 

    
%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/11/2022        1.9.1           hme
%
    
if nargin < 2
  tol=10^6*eps;
end
hd=harsanyi_dividends(v);
psQ=all(hd>=-tol);
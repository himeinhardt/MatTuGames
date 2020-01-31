function ccQ=AntiCoreCoverQ(v,tol)
% ANTICORECOVERQ checks if the anti-core cover a TU game v is non-empty.
%
% Usage:  ccQ=AntiCoreCoverQ(v,tol
% Define variables:
%  output:
%  ccQ      -- Returns 1 (true) whenever the core cover exists, 
%              otherwise 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- A tolerance value. Default is 10^7*eps
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/18/2015        0.7             hme
%                


if nargin < 2
  tol=10^7*eps;
end

ccQ=compromiseAntiAdmissibleQ(v,tol);

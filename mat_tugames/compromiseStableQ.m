function csQ=compromiseStableQ(v,tol)
% COMPROMISESTABLEQ checks if the game is compromise stable, that is,
% the core cover and the core coincide.
%
% Usage: csQ=compromiseStableQ(v,tol)
% Define variables:
%  output:
%  csQ      .. Returns 1 (true) whenever the core cover coincide 
%              with the core, otherwise 0 (false).
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
%   07/26/2014        0.5             hme
%                


if nargin < 2
  tol=10^7*eps;
end

[~,csQ]=compromiseAdmissibleQ(v,tol);


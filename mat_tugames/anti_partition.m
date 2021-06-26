function ap=anti_partition(ptn,n);
% ANTI_PARTITION computes from a partition its anti partition.
%
% Usage: ap=anti_partition(ptn,n);
%
% Define variables:
%  output:
%  ap       -- The anti partition of a partition of the player set N.
%  input:
%  ptn      -- A partition of the player set N.
%  n        -- The number of players involved.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/04/2013        0.5             hme
%                

N=2^n-1;
ap=N-ptn;


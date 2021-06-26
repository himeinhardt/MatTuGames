function C=getSymCostMatrix(n,rbd,method)
% GETSYMCOSTMATRIX computes a pseudo random symmetric cost matrix from the cardinality of the player set 
% and a upper bound value to specify the range from which the random number are drawn.
%
% Usage: C=getSymCostMatrix(n,rbd)
% Define variables:
%  output:
%  C        -- A symmetric cost matrix.
%
%  input:
%  n        -- number of actors involved. Position one is the source by default.
%  rbd      -- Upper bound of the range from which the random values are drawn.
%              The range is specified by [0,rbd].
%  method   -- String to specify if random floating number or integers 
%              are drawn. Permissible strings are 'float' or 'integer'.
%              Default is integer.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/12/2021        1.9        hme
%


if nargin < 2
   rbd=30;
   method='integer';
elseif nargin < 3
   method='integer';	
end	

if strcmp('integer',method)
   M = tril(randi([0 rbd],n),-1);
elseif strcmp('float',method)
   M = tril(max(rbd.*randn(n,n) + rbd,2),-1);
else
   M = tril(randi([0 rbd],n),-1);	
end	
C = M + M' + zeros(n);

function pnQ=PrenuclQ(v,x,tol)
% PRENUCLQ checks whether the (pre-)imputation x is the pre-nucleolus of game v.
% This function requires Matlab's Optimization toolbox.
%
%

% Define variables:
%  output:
%  pnQ      -- It returns 1 (true) or 0 (false).
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) 
%  tol      -- A positive tolerance value. The default is set to
%              10^6*eps'.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/09/2013        0.5             hme
%   04/02/2015        0.7             hme
%                

if nargin < 3
   tol=10^6*eps;
end    

pnQ=balancedCollectionQ(v,x,tol);

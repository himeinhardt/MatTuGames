function sh=ShapleyValueM(clv)
% SHAPLEY_VALUEM computes the Shapley value by the marginal contributions 
% of TU-game v. For n>10 use the function ShapleyValue() instead.
%
% Usage: sh=ShapleyValueM(clv)
%
% Define variables:
%  output:
%  sh       -- The Shapley-value of a TU-game v.
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
%   10/27/2012        0.3             hme
%                

n=clv.tuplayers;
Mgc=AllMarginalContributions(clv);
sh=sum(Mgc)/factorial(n);

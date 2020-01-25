function sh=p_ShapleyValueM(clv);
% P_SHAPLEY_VALUEM computes the Shapley value by the marginal contributions 
% of TU-game v. For n>10 use the function ShapleyValue() instead.
% Using Matlab's PCT.
%
% Usage: sh=p_ShapleyValueM(clv)
%
% Define variables:
%  output:
%  sh       -- The Shapley-value of a TU-game v.
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
%   10/30/2012        0.3              hme
%                

n=clv.tuplayers;
Mgc=p_AllMarginalContributions(clv);
sh=sum(Mgc)/factorial(n);

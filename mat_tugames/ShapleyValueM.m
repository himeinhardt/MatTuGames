function sh=ShapleyValueM(v)
% SHAPLEYVALUEM computes the Shapley value by the marginal contributions 
% of TU-game v. For n>10 use the function ShapleyValue() instead.
%
% Usage: sh=ShapleyValueM(v)
% Define variables:
%  output:
%  sh       -- The Shapley-value of a TU-game v.
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/10/2010        0.1 beta        hme
%   10/27/2012        0.3             hme
%                



N=length(v);
[~, n]=log2(N);

Mgc=AllMarginalContributions(v);
sh=sum(Mgc)/factorial(n);

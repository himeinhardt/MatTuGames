function sh=p_ShapleyValueM(v);
% P_SHAPLEY_VALUEM computes the Shapley value by the marginal contributions 
% of TU-game v. For n>10 use the function ShapleyValue() instead.
% Using Matlab's PCT.
%
% Usage: sh=p_ShapleyValueM(v)
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
%   05/20/2011        0.1 alpha        hme
%   10/27/2012        0.3              hme
%                



N=length(v);
[~, n]=log2(N);

Mgc=p_AllMarginalContributions(v);
sh=sum(Mgc)/factorial(n);

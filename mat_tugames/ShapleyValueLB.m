function [sh, cfc]=ShapleyValueLB(v)
% SHAPLEYVALUELB computes the Shapley value of a TU-game v.
% from a linear basis.
% For n>12 this function needs some time to complete.
%
% Usage: [sh cfc]=ShapleyValueLB(v)
%
% Define variables:
%  output:
%  sh       -- The Shapley-value of a TU-game v.
%  cfc      -- Coefficients of the linear basis of a game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%
%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/21/2013        0.4             hme
%

N=length(v);
[~, n]=log2(N);

if N==1
  sh=v;return;
 else
end

cfc=coeff_linearbasis(v);
k=1:n;
sC=2.^(k-1);
sh=cfc(sC);


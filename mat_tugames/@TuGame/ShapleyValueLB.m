function [sh, cfc]=ShapleyValueLB(clv)
% SHAPLEYVALUELB computes the Shapley value of a TU-game v.
% from a linear basis.
% For n>12 this function needs some time to complete.
%
% Usage: [sh cfc]=ShapleyValueLB(clv)
%
% Define variables:
%  output:
%  sh       -- The Shapley-value of a TU-game v.
%  cfc      -- Coeffient of the linear basis of a game v.
%
%  input:
%  clv      -- TuGame class object.
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

N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;

if N==1
  sh=v;return;
 else
end

cfc=coeff_linearbasis(clv);
k=1:n;
sC=2.^(k-1);
sh=cfc(sC);


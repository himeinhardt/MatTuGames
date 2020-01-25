function x=StandardSolution(v)
% STANDARDSOLUTION computes the standard solution of a two person game v.
%
% Usage: x=StandardSolution(v) 
% Define variables:
%  output:
%  x        -- The standard solution of a two person game v.
%
%  input:
%  v        -- A Tu-Game v of length 3.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/26/2010        0.1 beta        hme
%                

k=1:2;
x=v(k)+(v(3)-v(1)-v(2))/2;


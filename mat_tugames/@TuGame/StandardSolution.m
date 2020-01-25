function x=StandardSolution(clv)
% STANDARDSOLUTION computes the standard solution of a two person game v.
%
% Usage: x=StandardSolution(clv)
% 
% Define variables:
%  output:
%  x        -- The standard solution of a two person game v.
%
%  input:
%  clv      -- TuGame class object. (tuvalues of length 3)

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%                
v=clv.tuvalues;

k=1:2;
x=v(k)+(v(3)-v(1)-v(2))/2;


function x=sm_PreKernel(clv)
% SM_PREKERNEL computes an element of the simplified pre-kernel of game v.
%
% Usage: [x, fmin]=clv.sm_PreKernel(tol)
%
% Define variables:
%  output:
%  x         -- An element of the simplified pre-kernel of game v.
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
%   11/16/2017        0.9             hme

v=clv.tuvalues;
dv=clv.dual_game();
av=(v+dv)/2;
x=PreKernel(av);

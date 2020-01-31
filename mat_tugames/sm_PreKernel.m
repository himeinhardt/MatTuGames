function x=sm_PreKernel(v)
% SM_PREKERNEL computes an element of the simplified modified pre-kernel of game v.
%
% Usage: [x, fmin]=sm_PreKernel(v,tol)
%
% Define variables:
%  output:
%  x         -- An element of the simplified pre-kernel of game v.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/16/2017        0.9             hme

dv=dual_game(v);
av=(v+dv)/2;
x=PreKernel(av);

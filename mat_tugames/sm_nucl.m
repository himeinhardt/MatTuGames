function x=sm_nucl(v,tol)
% SM_NUCL computes the simplified modified nucleolus of game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=sm_nucl(v,tol)
% Define variables:
%  output:
%  x         -- The simplified nucleolus of game v.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/16/2017        0.9             hme

if nargin<2
 tol=10^6*eps;
end

dv=dual_game(v);
av=(v+dv)/2;
x=nucl_llp(av,tol);

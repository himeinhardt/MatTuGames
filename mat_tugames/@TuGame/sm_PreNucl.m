function x=sm_PreNucl(clv,tol)
% SM_NUCL computes the simplified pre-nucleolus of game v using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=clv.sm_PreNucl(tol)
% Define variables:
%  output:
%  x         -- The simplified pre-nucleolus of game v.
%
%  input:
%  clv      -- TuGame class object.
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

v=clv.tuvalues;
dv=clv.dual_game();
av=(v+dv)/2;
x=PreNucl_llp(av,tol);

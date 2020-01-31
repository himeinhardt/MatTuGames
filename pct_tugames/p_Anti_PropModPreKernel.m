function x=p_Anti_PropModPreKernel(v,x)
% P_ANTI_PROPMODPREKERNEL computes from (v,x) a proper modified anti-pre-kernel element using Matlab's PCT.
%
%  Inspiered by P. Sudh¨olter. Nonlinear Self Dual Solutions for TU-Games. In Potters J.A.M. Raghavan T.E.S. Ray D. Sen A.
%               Parthasarathy T., Dutta B., editor, Game Theoretical Applications to Economics and Operations Research, volume
%               18 of Theory and Decision Library: Series C, pages 33–50, Boston, MA, 1997b. Springer.
%
% Usage: [x Lerr smat xarr]=p_Anti_PropModPreKernel(v,x)
%
% Define variables:
%  output:
%  x        -- A proper modified anti-pre-Kernel element (output)
%  Lerr     -- List of computed function values of hx and h. 
%  smat     -- Matrix of maximum surpluses.
%  xarr     -- History of computed solution at each iteration step.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/10/2018        1.0             hme
%


N=length(v);
[~, n]=log2(N);

vm=DualFloor(v);
y=p_Anti_PreKernel(vm);
x=y(1:n);

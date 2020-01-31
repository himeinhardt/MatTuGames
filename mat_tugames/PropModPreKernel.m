function x=PropModPreKernel(v,x)
%PROPMODPREKERNEL computes from (v,x) a proper modified pre-kernel element.
%
%  Source: P. Sudh¨olter. Nonlinear Self Dual Solutions for TU-Games. In Potters J.A.M. Raghavan T.E.S. Ray D. Sen A.
%          Parthasarathy T., Dutta B., editor, Game Theoretical Applications to Economics and Operations Research, volume
%          18 of Theory and Decision Library: Series C, pages 33–50, Boston, MA, 1997b. Springer.
%
%          H. I. Meinhardt. Reconsidering Related Solutions of the Modiclus. Technical report, Karlsruhe Institute of Technology (KIT),
%          Karlsruhe, Germany, 2018. URL http://dx.doi.org/10.13140/RG.2.2.27739.82729.
%
%
%
% Usage: [x Lerr smat xarr]=PropModPreKernel(v,x)
%
% Define variables:
%  output:
%  x        -- A proper modified Pre-Kernel element (output)
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

vm=DualCover(v);
if nargin<2
  y=dc_PreKernel(vm);
else
  y=[x,x];
  y=dc_PreKernel(vm,y); 
end
x=y(1:n);

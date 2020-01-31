function MPKQ=PrkEqsModPrkQ(v,x,str,tol)
% PRKEQSMODPRKQ checks whether a pre.kernel element is also an element of the modified
% as well as proper modified pre-kernel. 
%
%  Source: H. I. Meinhardt. Reconsidering Related Solutions of the Modiclus. Technical report, Karlsruhe Institute of Technology (KIT),
%          Karlsruhe, Germany, 2018c. URL http://dx.doi.org/10.13140/RG.2.2.27739.82729. 
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          P. Sudhoelter. Nonlinear Self Dual Solutions for TU-Games. In Potters J.A.M. Raghavan T.E.S. Ray D. Sen A.
%          Parthasarathy T., Dutta B., editor, Game Theoretical Applications to Economics and Operations Research, volume
%          18 of Theory and Decision Library: Series C, pages 33–50, Boston, MA, 1997b. Springer. 
% 
%  Usage: MPKQ=PrkEqsModPrkQ(v,x,'ecc')
%
%
% Define variables:
%
%  output: Fields
% mpkeqQ       -- Returns 1 (true) whenever x is a pre-kernel as well as a 
%                 modified pre-kernel element, otherwise 0 (false).
% pmpkeqQ      -- Returns 1 (true) whenever x is a pre-kernel as well as a 
%                 proper modified pre-kernel element, otherwise 0 (false).
% smat         -- Returns the matrix of maximum surpluses.
% smat_mpk     -- Returns the matrix of maximum modified surpluses.
% smat_pmpk    -- Returns the matrix of maximum proper modified surpluses.
% led          -- Returns 1 whenever, the game fulfills LED w.r.t. x,
%                 otherwise 0 (false).
% lex          -- largest excess.
% ldex         -- largest dual excess.
% ex           -- excess array.
% dex          -- dual excess array.
% vx           -- Original or ECC depending on the game context 
%                 provided by str.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient. 
%              Default is the pre-kernel.
%  str      -- A string that defines different Methods.
%              Permissible methods are:
%              'orig' that is, game v is referred to.
%              'ecc' that is, excess comparability cover game vx is referred to. 
%              Default is 'ecc'.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional)
%    

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/23/2014        1.0             hme
%


if nargin<1
  error('A game must be given!')
elseif nargin<2
 x=PreKernel(v);  
 tol=10^6*eps;
else
  tol=10^6*eps;   
end

N=length(v);
[~, n]=log2(N);

if strcmp(str,'orig')
   vx=v;
elseif strcmp(str,'ecc')
   vx=ECCoverGame(v,x);
else
   vx=ECCoverGame(v,x);
end


dvx=dual_game(vx);
LED=LED_propertyQ(vx,x);
[pkQ, smat]=PrekernelQ(vx,x);
[mpkQ, smat_mpk]=ModPrekernelQ(vx,x);
[pmpkQ, smat_pmpk]=PModPrekernelQ(vx,x);

mpkglQ=all(all(abs(smat-smat_mpk)<=tol));
mpkeqQ=all([pkQ,mpkQ,mpkglQ]);
ex=excess(vx,x);
dex=excess(dvx,x);
lex=max(ex);
ldex=max(dex);
mt1=ones(n)*max(smat(:));
Ar=[smat,mt1;mt1,smat];
%%abs(Ar-smat_pmpk)<=tol
pmpkglQ=all(all(abs(Ar-smat_pmpk)<=tol));
pmpkeqQ=all([pkQ,pmpkQ,mpkglQ]);
%% Selecting only the relevant range for smat_pmpk
midx=[true(n,n),false(n,n);false(n,n),false(n,n)];
smat_pmpk_rs=reshape(smat_pmpk(midx),n,n);

MPKQ.mpkeqQ=mpkeqQ;
MPKQ.pmpkeqQ=pmpkeqQ;
MPKQ.smat=smat;
MPKQ.smat_mpk=smat_mpk;
MPKQ.smat_pmpk=smat_pmpk_rs;
MPKQ.led=LED;
MPKQ.lex=lex;
MPKQ.ldex=ldex;
MPKQ.ex=ex;
MPKQ.dex=dex;
MPKQ.vx=vx;


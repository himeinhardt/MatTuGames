function PD=GetProbDist(v,x,tol)
% GETPROBDIST tries to find for the game v the corresponding
% probability distribution over the players' orderings w.r.t. pre-imputation x.  
% 
%  Usage: PRKS=GetProbDist(v,x)
%
%  Inspired by R. J. Weber (1988), Probabilistic values for game, In The Shapley Value. Essays in Honor of Lloyd Shapley. 
%                 Pierre Dehez (2017), On Harsanyi Dividends and Asymmetric Values. IGTR, Vol. 19. No. 3, 36 pages. 
%
%
% Define variables:
%  output: Fields
%  pS       -- Returns the probability distribution over the players' orderings.
%              Shapley value must return the uniform probability distribution.
%  prbQ     -- Returns true, if pS is a probability distribution, otherwise zero.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- Payoff vector of size(1,n), must be an element of the Weber set.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/12/2022        1.9.1           hme
%


if nargin<2
   x=PreKernel(v);
   tol=10^6*eps;
elseif nargin < 3
   tol=10^6*eps;
end	

mrg=AllMarginalContributions(v)';
x=x';
pS=pinv(mrg)*x;
prb1Q=abs(sum(pS)-1)<=tol;
prb2Q=all(pS>=-tol);
PD.pS=pS';
PD.prbQ=prb1Q & prb2Q;

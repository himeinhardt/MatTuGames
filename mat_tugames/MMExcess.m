function MMEX=MMExcess(v,x)
% EXCESS computes the minimal and maximal excess vector of game v 
% and its dual w.r.t. x.
%
% Usage: mmex=MMExcess(v,x)
%
% Define structure variable:
%
% output:
% mev        -- maximal excess vector of game v w.r.t. x.
% medv       -- maximal excess vector of the dual game of v w.r.t. x.
% sev        -- minimal excess vector of game v w.r.t. x.
% sedv       -- minimal excess vector of the dual game of v w.r.t. x.
% x          -- payoff (input vector) of size (1,n).
%
% input: 
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/19/2018        1.0             hme


ex=excess(v,x);
dv=dual_game(v);
dex=excess(dv,x);

MMEX.mev=max(ex);
MMEX.medv=max(dex);
MMEX.sev=min(ex);
MMEX.sedv=min(dex);
MMEX.x=x;

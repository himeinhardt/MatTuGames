function [v,rv]=p_GetMarketGame(n,seed,tol)
% P_GETMARKETGAME determines from a random generated game the corresponding market game
% using MATLAB's PCT.
%
% Source: Shapley and Shubik, On Market Games, JET 1, 9-25 (1969).
%
% Usage: v=p_GetMarketGame(v,'seed'tol)
%
% Define output variable:
%  v        -- A Market game of size 2^n-1.
%  rv       -- Random Game of size 2^n-1. 
%
%  input:
%  n        -- Size of the player set, must be an integer. Default value is n=4.
%  seed     -- A seed value to randomly generate a TU game. Default value is seed=137.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/09/2022        1.9.1           hme
%


if nargin < 1
   n=4;
   seed=137;
   tol=10^8*eps;
elseif nargin < 2
   seed=137;
   tol=10^8*eps;
elseif nargin < 3
   tol=10^8*eps;  	
end	
N=2^n-1;
rng(seed,'v5uniform');
rv=ceil(100*rand(1,N));
[~,w]=p_totallyBalancedCoverQ(rv);
[tcbQ,w]=p_totallyBalancedCoverQ(w,tol);
if tcbQ==1
   v=w;
else
   warning('GMarket:Exit','No Market Game determined! Change the numerical tolerance value.')	
   v=[];  	
end	


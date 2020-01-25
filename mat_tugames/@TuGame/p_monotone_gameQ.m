function [mq A]=p_monotone_gameQ(clv,tol)
% P_MONOTONE_GAMEQ returns 1 whenever the game v is monotone.
% Using Matlab's PCT.
%
% Usage: [mq A]=p_monotone_gameQ(clv,tol)
%
% Define variables:
%  output:
%  mq       -- Returns 1 (true) or 0 (false).
%  A        -- A vector which indicates whether the 
%              marginal contribution of player i are 
%              increasing (1) or not (0).
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. By default, it is set to (-2*10^4*eps).
%              (optional) 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/30/2012        0.3              hme
%                

if nargin<2
   tol=(-2*10^4*eps);
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

A=false(n,1);
S=1:N;

parfor i=1:n
   zi=bitset(0,i); % empty set.
   dz=v(zi)>=tol; 
   Tni=bitget(S,i)==0; % non-empty sets. 
   T=S(Tni);
   Ti=bitset(T,i);
   dvi=v(Ti)-v(T)>=tol;
   A(i)=all([all(dvi),dz]);
end

mq=all(A);

function [acvQ acvq]=p_AlmostConvex_gameQ(clv,tol)
% P_ALMOSTCONVEX_GAMEQ returns 1 whenever the game v is almost convex using Matlab's PCT.
%
% Source: The bargaining set and the kernel for almost-convex games
%
% Usage: [acvQ acvq]=clv.p_AlmostConvex_gameQ(tol)
%
% Define variables:
%  output:
%  acvQ     -- Returns 1 (true) whenever all proper sub-games are convex, otherwise 0 (false).
%  acvq     -- Returns an array of 1/0 to indicate whether the sub-game S is convex. 
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
%   06/05/2020        1.9             hme
%                

if nargin<2
   tol=-2*10^4*eps;
end

N=clv.tusize;
n=clv.tuplayers;
N1=N-1;
acvQ=false(1,N1);

parfor S=1:N1
    sg=clv.SubGame(S);
    acvq(S)=convex_gameQ(sg,tol);
end	
acvQ=all(acvq);


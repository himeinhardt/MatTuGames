function adQ=AdditiveQ(clv,tol)
% ADDITIVEQ checks if the game v is additive or not. Notice that some
% authors denote this propertiy as essential.
%
% Usage: eQ=EssentialQ(v) 
%
% Define variables:
%  output:     
%  adQ       -- Returns 1 (ture) whenever the game is additive,
%              otherwise 0 (false)
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- A tolerance value, defualt is tol=10^6*eps.
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/21/2014        0.5             hme

if nargin <2
  tol=10^6*eps;
end

N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;

k=1:n;
vk=v(bitset(0,k));
av=additive_game(vk);
adQ=all(abs(av-v)<tol);

function ex=excess(clv,x)
% EXCESS computes the excess vector of game v w.r.t. x.
%
% Usage: ex=clv.excess(x)
%
% Define variables:
% output:
% ex         -- excess vector of game v w.r.t. x.
%
% input: 
%  clv        -- TuGame class object.
%  x        -- payoff vector of size(1,n).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%                

% Computing the excess vector w.r.t. x.
v=clv.tuvalues;
n=clv.tuplayers;
% Borrowed from J. Derks
Xm=x(1); for ii=2:n, Xm=[Xm x(ii) Xm+x(ii)]; end
ex=v-Xm;

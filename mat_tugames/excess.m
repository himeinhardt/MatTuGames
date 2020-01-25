function ex=excess(v,x)
% EXCESS computes the excess vector of game v w.r.t. x.
%
% Usage: ex=excess(v,x)
% Define variables:
% output:
% ex         -- excess vector of game v w.r.t. x.
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
%   08/19/2010        0.1 beta        hme
%   06/08/2012        0.2 beta        hme
%   09/12/2012        0.2             hme
%                

n=length(x);
% Computing the excess vector w.r.t. x.
% Borrowed from J. Derks
Xm=x(1); for ii=2:n, Xm=[Xm x(ii) Xm+x(ii)]; end
ex=v-Xm;

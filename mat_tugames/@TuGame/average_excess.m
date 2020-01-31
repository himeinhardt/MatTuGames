function avex=average_excess(clv)
% AVERAGE_EXCESS computes the average excess of game v.
%
% Usage: avex=clv.average_excess()
% Define variables:
% output:
% avex       -- average excess of game v.
%
% input: 
%  clv        -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/22/2017        0.9             hme
% 

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
avex=(sum(v)-2^(n-1)*v(N))/N;

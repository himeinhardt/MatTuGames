function avex=average_excess(v)
% AVERAGE_EXCESS computes the average excess of game v.
%
% Usage: avex=average_excess(v)
% Define variables:
% output:
% avex       -- average excess of game v.
%
% input: 
%  v        -- A Tu-Game v of length 2^n-1. 
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


N=length(v);
[~, n]=log2(N);
avex=(sum(v)-2^(n-1)*v(N))/N;

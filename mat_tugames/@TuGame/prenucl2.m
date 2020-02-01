function [x,qq,pp,y]=prenucl2(clv,x)
%PRENUCL computes the prenucleolus of a game
% [x,qq,pp]=prenucl(v) computes the prenucleolus x of the game v.
% the output qq is the number of LP's, and pp the number of pivots.

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/27/04$
%   Jean Derks

N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;
F=N;
qq=1;
global TOLR;
if isempty(TOLR), TOLR=10^5*eps;end
%TOLR
[x,U,pp,t,y]=leastcore(clv,x);
if U==-1, return; end % no independent feasible coalitions
while 1  
    F=Ionly(F,U); % F stays 'independent' in this way?
    d=size(F);
    if d(1,2)==n, break; end
    qq=qq+1;
    [x,U,p,t,y]=leastcore(clv,x,F);
    if U==-1, return; end % no independent feasible coalitions
    pp=pp+p;
end

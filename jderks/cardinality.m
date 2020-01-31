function [n,r]=cardinality(x)
%CARDINALITY returns the number of players involved
% n=cardinality(x) n is the cardinality of x in case x is a number which
% represents a coalition. In case x is a matrix it represents a game, and
% in that case the (minimal) number of players is returned in this game.
% r returns the type of x: 0,1,2 if restricted, non-restr.,coalition

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/24/02$
%   Jean Derks

d=size(x);
if d(1,2)==1   % coalition
    S=x(1,1);
    l=lplayer(S);
    n=0;
    for i=1:l
    if bitget(S,i), n=n+1; end
    end 
    r=2;   
elseif d(1,1)>1 % restricted
    n=x(2,1);   
    r=0; 
else            % full game
    S=d(1,2)-1;
    n=lplayer(S);
    r=1;  
end;    

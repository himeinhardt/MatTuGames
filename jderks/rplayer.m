function n=rplayer(S)
%RPLAYER index of the rightmost player in S
% n=rplayer(S) computes the smallest index n of a player
% in S 

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/24/02$
%   Jean Derks

n=1; 
while bitget(S,n)==0, n=n+1; end

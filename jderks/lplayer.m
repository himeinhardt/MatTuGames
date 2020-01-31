function n=lplayer(S)
%LPLAYER index of the leftmost player in S
% n=lplayer(S) computes the smallest number n of players
% such that S is a coalition in an n-person game

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/24/02$
%   Jean Derks

n=0; 
while S>bitshift(1,n), n=n+1; end
if S==bitshift(1,n), n=n+1; end;

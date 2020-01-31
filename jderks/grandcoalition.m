function N=grandcoalition(n)
%GRANDCOALITION returns a coalition of max cardinality
% N=grandcoalition(n) N is a coalition of cardinality n, if n single;
% otherwise n is considered a game and N will hold its grand coalition

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.01$  $Date: 2005/24/02$
%   Jean Derks

d=size(n);
k=d(1,2);
if k==1, k=n(1,1); 
else
    k=lplayer(k);
end
N=2^k-1;
%N=bitcmp(0,k)

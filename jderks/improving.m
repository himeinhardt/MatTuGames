function [y,i]=improving(F,U,n)
%IMPROVING returns an allocation, gaining 1 for some coalitions, zero for others.
% y=mproving(F,U,n) allocation y is computed such that the coalitions in F gain zero, and y(S)=1
% for the coalitions S in U.

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/27/02$
%   Jean Derks

global TOLR;
if isempty(TOLR), TOLR=100*eps; end

d=size(F);
mf=d(1,2);
d=size(U);
mu=d(1,2);
m=mf+mu;
G=zeros(m+1,n+1);
for i=1:mu, for j=1:n, G(i,j)=bitget(U(i),j); end, end
j=n+1;
for i=1:mu, G(i,j)=1; end
for i=1:mf, for j=1:n, G(mu+i,j)=bitget(F(i),j); end, end

G(m+1,n+1)=1;
y=UT(G',1);

if y(1)==Inf % case 3
    i=0;
else
    
    for i=1:mu, if y(i)<=-TOLR, break; end, end
    if y(i)<=-TOLR,  for j=1:n+1, G(i,j)=0; end   % case 2
    else
        i=-1; % case 1
        return;
    end
end
G(m+1,n+1)=0;
y=UT(G,1)';


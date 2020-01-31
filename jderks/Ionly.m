function F=Ionly(F1,F2)
%IONLY selects the independent coalitions in a collection
% F=Ionly(F1,F2) selects the independent 
% coalitions in the collection F1, or the union of the two collections F1 and F2.

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/27/04$
%   Jean Derks


d=size(F1);
mf=d(1,2);
S=F1(1);
for i=2:mf, S=bitor(F1(i),S); end
if nargin==2
    d=size(F2);
    mu=d(1,2);
    for i=1:mu, S=bitor(F2(i),S); end
else
    mu=0;
end
n=lplayer(S);
m=mf+mu;
G=zeros(m,n);
for i=1:mu, for j=1:n, G(i,j)=bitget(F2(i),j); end, end
for i=1:mf, for j=1:n, G(mu+i,j)=bitget(F1(i),j); end, end

F=UT(G,0); % supplies the indices of the rows that constitute a max. indep. set 

d=size(F);
m=d(1,2);
for i=1:m 
    if F(i)>mu
        F(i)=F1(F(i)-mu);
    else
        F(i)=F2(F(i));
    end
end


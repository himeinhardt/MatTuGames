function U=UT(H,t)
%UT computes a squared triangulated version of a matrix
% U=UT(H,t) computes a squared triangulated version of H.
% if t is supplied, it determines the output:
%       t==0: U is the vector of row indices, together determining the
%               independency;
%       t>0: the output is a kernel element of the column space, with -1 on 
%               the last position (this coordinate is not supplied). An example:
%               the expression x=UT([A b],n+1) provides a solution of Ax=b 
%               (x=Inf in case unsolvable))

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/22/03$
%   Jean Derks

global TOLR;
if isempty(TOLR), TOLR=100*eps; end

G=H;
d=size(G);
n=d(1,2);
m=d(1,1);
M=zeros(1,m);
U=zeros(n,n);

for j=1:n   % search for a nonzero entry in column j
    for i=1:m, if G(i,j)>TOLR | G(i,j)<-TOLR, break; end, end
    if G(i,j)>TOLR | G(i,j)<-TOLR
        M(i)=1;
        p=1/G(i,j);
        U(j,j)=1;
        for k=j+1:n
            U(j,k)=p*G(i,k); 
            G(i,k)=0;
        end
    end
    for ii=i+1:m
        if G(ii,j)>TOLR | G(ii,j)<-TOLR
            p=G(ii,j);
            for k=j+1:n, G(ii,k)=G(ii,k)-U(j,k)*p; end
        end
    end
end

if nargin==1, return; end
if t==0
    u=sum(M);
    U=zeros(1,u);
    i=1;
    while u>0
        if M(i) 
            U(u)=i;
            u=u-1;
        end
        i=i+1;
    end
    return;
end

if U(n,n)>TOLR | U(n,n)<-TOLR
       U=Inf;
       return;
end

for i=1:n-1 
    if (U(i,i)<=TOLR & U(i,i)>=-TOLR) & (U(i,n)>TOLR | U(i,n)<-TOLR)
       U=Inf;
       return;
    end
end
   
    
M=zeros(n-1,1);
for i=n-1:-1:1, 
    if (U(i,i)>TOLR | U(i,i)<-TOLR)
        t=U(i,n);
        for j=i+1:n-1, t=t-M(j)*U(i,j); end
        M(i)=t;    
    end
end
U=M;

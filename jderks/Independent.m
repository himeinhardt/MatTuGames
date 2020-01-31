function [IND,U]=Independent(F,n,ind)
%INDEPENDENT computes the data of independent coalitions 
% [IND,U]=Independent(F,n) if F is a set of coalitions on a player set of 
% cardinality n, then the function computes for each coalition whetheer its 
% indicator function can be written as a combination of the indicator functions 
% of the coalitions in F. These 2^n-1 outcomes are stored in the array IND,as a sequence
% of 0 and 1's. See compact() and expand() for details on this storage issue.
% The output U is the coalition of players i, such that adding the coalitions {i} to F
% evolves a collection of full rank (cardinality of U is therefore n-rank(F)).

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/22/03$
%   Jean Derks

global TOLR;
if isempty(TOLR), TOLR=100*eps; end

d=size(F);
f=d(1,2);
U=F(1);
for i=2:f, U=bitor(U,F(i)); end
if nargin<2, n=lplayer(U); end
N=grandcoalition(n);

G=zeros(f,n);
for i=1:f, for j=1:n, if bitget(F(i),j), G(i,n-j+1)=1; end, end, end
G=UT(G); % G is now a square matrix, extracted from F by triangulation

U=0;
for i=1:n
    if G(i,i)==0
        U=bitset(U,i);
        G(i,i)=1;
    end
end

Inc=zeros(n,n);
for i=1:n
    Inc(i,n-i+1)=1;
    for j=n-i+2:n, Inc(i,j)=-1; end
end

% cnt=0; test of obtaining all vectors e^S in about 2^(n+1) steps
% s=zeros(1,n);
% s(n)=1;
% for S=1:N-1
%     k=rplayer(N-S);
%     for j=n-k+1:n
%        cnt=cnt+1;
%        s(j)=s(j)+Inc(k,j);
%    end
% end
% cnt

if nargin>2 % restricted game
    d=size(ind); % ind() feasible coalitions
    m=d(1,2);
    IND=zeros(1,m);
    W=inv(G);
    plu=players(U);
    u=cardinality(U);
    for i=1:m
        S=ind(i);
        if cardinality(S)==1
            lambda=W(n+1-players(S),:);
        else
            lambda=sum(W(n+1-players(S),:));
        end
        for j=1:u
            if lambda(plu(j))<-TOLR | lambda(plu(j))>TOLR 
                IND(i)=1; 
                break; 
            end
        end
    end
    return;
end

W=Inc*inv(G); % is again an upper triangular m. like Inc

% now checking all independencies, and stores the info in a compactified simple game IND.
% See compact()

p=52;
P=floor(N/52);
IND=zeros(1,P+1);

lambda=W(1,:); 
if bitget(U,n)
    t=n; % leftmost index in lambda with nonzero entry and bitget(U,k)
    IND(1)=1; % independency
else 
    t=n+1; 
end

Q=1;q=1; % indices for storing independency in the data IND
for S=1:N-1
    q=q+1; 
    if q>p
        q=1;
        Q=Q+1;
    end
    
     k=rplayer(N-S);
     for j=n-k+1:n
        lambda(j)=lambda(j)+W(k,j);
    end
    % lambda*G is the indicator function of S+1
    % If there are nonzero entries in lambda wrt. indices in U
    % (the added players to F, for obtaining the non-singular G)
    % then S is independent of the coalitions in F
    % We cannot just test all entries since it destroys the speed factor of 2^(n+1)
    if t<n-k+1 
        IND(Q)=bitset(IND(Q),q); % independency
    else
        for j=n-k+1:n
            if (lambda(j)>TOLR | lambda(j)<-TOLR) & bitget(U,j)~=0
                t=0; 
                break;
            end
        end
        if t==0
            t=j;
            IND(Q)=bitset(IND(Q),q); % independency
        else
            t=n+1;
        end
    end
end

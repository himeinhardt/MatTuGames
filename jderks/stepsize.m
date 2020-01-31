function [lambda,U]=stepsize(x,y,t,v,IND)
%STEPSIZE computes the optimal stepsize in a pivot of the least core algorithm

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/27/04$
%   Jean Derks

global TOLR;
if isempty(TOLR), TOLR=100*eps; end

n=cardinality(v);
N=grandcoalition(n);
upto_1=1-TOLR;

d=size(v);
if d(1,1)>1 % restricted
    m=d(1,2);
    U=0;
    lambda=flintmax;
    for i=2:m;
       S=v(2,i);
       plS=players(S);
       sy=sum(y(plS));
       if IND(i) & sy<upto_1
           lda=(sum(x(plS))+t-v(1,i))/(1-sy);
           if lda<lambda
               U=S;
               lambda=lda;
           end
       end
    end
    return;
end


incx=zeros(1,n);
incy=zeros(1,n);

for i=1:n
    lambda=x(i);
    for j=1:i-1, lambda=lambda-x(j); end
    incx(i)=lambda;
    lambda=y(i);
    for j=1:i-1, lambda=lambda-y(j); end
    incy(i)=lambda;
end

sx=incx(1); % testing with x=[.00001,.0001,.001,.01,.1]
sy=incy(1);  % y=[1          10         100        1000       10000]
Q=1;
q=1;
p=52;
if bitget(IND(Q),q) & sy<upto_1
    U=1;
    lambda=(sx+t-v(1))/(1-sy);
else
    U=0;
    lambda=flintmax;
end
for S=1:N-1

    q=q+1;
    if q>p
        q=1;
        Q=Q+1;
    end

    k=rplayer(N-S);
    sx=sx+incx(k);
    sy=sy+incy(k);

    if bitget(IND(Q),q) & sy<upto_1
        lda=(sx+t-v(S+1))/(1-sy);
        if lda<lambda
            U=S+1;
            lambda=lda;
        end
    end
end


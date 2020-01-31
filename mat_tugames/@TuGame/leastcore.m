function [x,U,p,t,y]=leastcore(clv,xb,F,G)
%LEASTCORE computes an allocation of the least core
% [x,U,p,t]=leastcore(v,xb,F,G) for a game v the outcome will be 
% an allocation x from the least core. If supplied, the argument xb
% will serve as a starting point (efficiency is assumed), and if 
% a collection of (independent) coalitions F is supplied, only values are considered
% of coalitions that are independent of those in F. And if furthermore G is supplied
% then x should obey x(S)\geq v(S) for all coalitions S in the collection G.

% The extra output U is a set of coalitions S, independent of those in F, and 
% for which all allocations x in the least core attain the same outcome, x(S)=v(S)+t.
% The output p is the number of pivots.
% And t is the value, with which each coalition value should be adapted in order to 
% attain the game with core equal to the least core of v.

%   Copyright 2005 Universiteit Maastricht, dept. of Mathematics
%   $Revision: 1.00$  $Date: 2005/27/04$
%   Jean Derks

global TOLR;
if isempty(TOLR), TOLR=100*eps; end 

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

if nargin<2
    gvN=Rval(v,N)/n;   
    x=ones(1,n)*gvN; 
else
    x=xb;
end
if nargin<3, F=[N]; end
if nargin==4 
    w=shrink(v,G);      % the G-restriction of v
    if unstable(v,x,G)>0
        x=leastcore(w);
        if unstable(w,x)    % the G-restricted core of v is empty
            x=-1;           % no solution possible
            return
        end
    end
end

d=size(v);
if d(1,1)>1 % restricted game
    IND=Independent(F,n,v(2,:));
    if sum(IND)==0 % all feasible coalitions are dependent; so, each efficient solution is ok.
        U=-1;
        return;
    end 
else
    IND=Independent(F,n);
end

p=0;
[lambda,U]=stepsize(x,zeros(1,n),0,v,IND);
if U==0, return; end % no improvement possible

t=-lambda;
while 1
    p=p+1;
    [y,i]=improving(F,U,n); % y zero for the coalitions in F (yielding y(N)=0), and y(U)=1
    
    if i<0
        d=size(U);
        for k=1:d(1,2), if y(k)<=TOLR, U(k)=0; end, end  % tolr is needed, in stead of just 0...
        break; 
    end  % no improving direction available, thus x optimal

    if y==Inf % error!!
        if d(1,1)>1 
            h=-1
        else
            warning('PRN:Error','Error triggers! Probably no pre-nucleolus found! Check solution!')
            break;
            s=expand(IND);
            f=length(F);
            h=zeros(1,f);
            for k=1:f, h(k)=s(F(k)); end
            h
        end
        return;
    end
    
    [lambda,S]=stepsize(x,y,t,v,IND);
    if S==0, return; end % not to be expected

    x=x+lambda*y;
    t=t-lambda;
    if i>0, U(i)=S; % U(i) free spot (see improving...)
    else
        U=[U, S];
    end  
end

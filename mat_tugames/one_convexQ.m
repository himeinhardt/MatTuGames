function ocvQ=one_convexQ(v,tol)
% ONE_CONVEXQ checks whether the game v is 1-convex.
%
%
% Usage: ocvQ=one_convexQ(v)
% Define variables:
%  output:
%  ocvQ     -- Returns a one if the game is 1-convex, otherwise zero.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- A tolerance value. Default is set to tol=10^6*eps;
%

%    
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/16/2022        1.9.1           hme
%

if nargin<2    
   tol=10^6*eps; 
end    
N=length(v);
[~, n]=log2(N);
bv=zeros(1,n);


[~,bv]=Gap(v);    

grQ=sum(bv) + tol >=v(N);
pl=1:n;
for k=1:N-1
    aS=bitget(k,pl)==0;
    grSQ(k)=v(k)+sum(bv(aS)) <=v(N) +tol; 
end   
agrSQ=all(grSQ);
ocvQ= grQ & agrSQ;

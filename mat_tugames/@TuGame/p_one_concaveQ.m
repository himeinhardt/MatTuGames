function ocvQ=p_one_concaveQ(clv,tol)
% P_ONE_CONCAVEQ checks whether the game v is 1-concave while 
% using MATLAB's PCT.
%
%
% Usage: ocvQ=clv.p_one_concaveQ()
% Define variables:
%  output:
%  ocvQ     -- Returns a one if the game is 1-concave, otherwise zero.
%
%  input:
%  clv      -- TuGame class object.
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
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
bv=zeros(1,n);


[~,bv]=clv.p_Anti_Gap();    

grQ=sum(bv)<=v(N) + tol;
pl=1:n;
parfor k=1:N-1
    aS=bitget(k,pl)==0;
    grSQ(k)=v(k)+sum(bv(aS)) +tol>=v(N); 
end   
agrSQ=all(grSQ);
ocvQ= grQ & agrSQ;

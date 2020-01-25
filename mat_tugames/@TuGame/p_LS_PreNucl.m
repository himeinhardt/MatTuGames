function x=p_LS_PreNucl(clv)
% P_LS_PRENUCL computes the least square pre-nucleolus of a game
% using Matlab's PCT.
%
% Usage: x=p_LS_PreNucl(clv)
% Define variables:
%  output:
%  x        -- Least square pre-nucleolus 
%
%  input:
%  clv      -- TuGame class object.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/10/2013        0.4             hme
%     
    

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
    
S=1:N;
sV=zeros(1,n);

 parfor k=1:n
     a=S(bitget(S,k)==1);
     sV(k)=sum(v(a));
 end     
nV=sum(sV);
sN=2^(n-2);
x=(v(N)/n)+(n*sV-nV)/(n*sN);

function sdsh=p_SD_ShapleyValue(v)
% P_SD_SHAPLEY_VALUE computes the surplus division Shapley value that
% is identical to the Shapley value of a TU-game v using Matlab's PCT.
%
% Source: David Pérez-Castrillo and David Wettstein (2001), JET 
% Bidding for the Surplus : A Non-cooperative Approach to the Shapley Value
%
%
% Usage: sdsh=SD_ShapleyValue(v)
% Define variables:
%  output:
%  sh       -- The surplus division Shapley-value of a TU-game v.
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/06/2018        1.0             hme
%

N=length(v);
[~, n]=log2(N);
if N==1
  sh=v;return;
 else
end


k=1:n;
Nk=bitset(N,k,0);
shg=zeros(n,n);
parfor ii=1:n
    sol=zeros(1,n);
    subg=SubGame(v,Nk(ii));
    idx=k(k~=ii);
    sol(idx)=ShapleyValue(subg);
    shg(ii,:)=sol;
end    
vN=v(N)-v(Nk);    
sdsh=(vN+sum(shg,1))/n;


    

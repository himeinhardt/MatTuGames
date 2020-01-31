function sdsh=p_SD_ShapleyValue(clv)
% P_SD_SHAPLEY_VALUE computes the surplus division Shapley value that
% is identical to the Shapley value of a TU-game v using Matlab's PCT.
%
% Source: David Pérez-Castrillo and David Wettstein (2001), JET 
% Bidding for the Surplus : A Non-cooperative Approach to the Shapley Value
%
%
% Usage: sdsh=clv.SD_ShapleyValue()
% Define variables:
%  output:
%  sh       -- The surplus division Shapley-value of a TU-game v.
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
%   04/06/2018        1.0             hme
%

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
vn=v(N);
Nk=clv.tuSi;
vNk=v(Nk);
if N==1
  sh=v;return;
 else
end
clear v;


shg=zeros(n,n);
parfor ii=1:n
    sol=zeros(1,n);
    subg=clv.SubGame(Nk(ii));
    idx=k(k~=ii);
    sol(idx)=ShapleyValue(subg);
    shg(ii,:)=sol;
end    
vN=vn-vNk;    
sdsh=(vN+sum(shg,1))/n;


    

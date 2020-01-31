function sdsh=SD_ShapleyValue(clv)
% SD_SHAPLEY_VALUE computes the surplus division Shapley value that
% is identical to the Shapley value of a TU-game v.
%
% Source: 
%
%
% Usage: sdsh=clv.SD_ShapleyValue()
% Define variables:
%  output:
%  sh       -- The Shapley-value of a TU-game v.
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
for ii=1:n
    subg=clv.SubGame(Nk(ii));
    idx=k(k~=ii);
    shg(ii,idx)=ShapleyValue(subg);
end    
vN=vn-vNk;    
sdsh=(vN+sum(shg,1))/n;


    

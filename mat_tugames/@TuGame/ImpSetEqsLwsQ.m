function ilQ=ImpSetEqsLwsQ(clv)
% IMPSETEQSQ checks if the lower set and the imputation set are equal.
%
%  Usage: ilQ=ImpSetEqsLwsQ(v) 
%
%
% Define variables:
%  output:
%  ilQ      -- Returns true (1) if the lower set and the imputation set coincide,
%              otherwise false (0).
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
%   12/23/2014        0.6             hme
%


v=clv.tuvalues;
N=clv.tusize;

tol=10^6*eps;
[impv,~,imp_vol,~]=clv.CddImputationVertices;
y1=range(impv);
[~,idx]=min(y1);
sma=clv.smallest_amount;
nlwseQ=sum(sma)>v(N);
if nlwseQ==1
  lws_vol=-Inf;
else
  [~,~,lws_vol,~]=clv.CddLowerSetVertices(idx);
end
ilQ=abs(lws_vol-imp_vol)<tol;

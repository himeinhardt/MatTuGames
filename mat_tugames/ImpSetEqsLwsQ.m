function ilQ=ImpSetEqsLwsQ(v)
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
%  v        -- A Tu-Game v of length 2^n-1. 
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


N=length(v);
[~, n]=log2(N);

tol=10^6*eps;
[impv,~,imp_vol,~]=CddImputationVertices(v);
y1=range(impv);
[~,idx]=min(y1);
sma=smallest_amount(v);
nlwseQ=sum(sma)>v(N);
if nlwseQ==1
  lws_vol=-Inf;
else
  [~,~,lws_vol,~]=CddLowerSetVertices(v,idx);
end
ilQ=abs(lws_vol-imp_vol)<tol;

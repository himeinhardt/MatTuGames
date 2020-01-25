function ruQ=ReasSetEqsUpsQ(v)
% REASSETEQSUPS checks if the reasonable set and the upper set are equal.
%
%  Usage: ruQ=ReasSetEqsUpsQ(v) 
%
%
% Define variables:
%  output:
%  ilQ      -- Returns true (1) if the upper set and the reasonable set coincide,
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
[reasv,~,reas_vol,~]=CddReasonableSetVertices(v);
y1=range(reasv);
[~,idx]=min(y1);
pa=proper_amount(v);
usQ=sum(pa)>=v(N);
if usQ==1
  [~,~,ups_vol,~]=CddUpperSetVertices(v,idx);
else
  ups_vol=-Inf;
end
ruQ=abs(ups_vol-reas_vol)<tol;

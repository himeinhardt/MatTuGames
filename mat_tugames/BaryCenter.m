function cnt=BaryCenter(v)
% BaryCENTER computes the barycenter of the core.
%
% Usage: cnt=BaryCenter(v)
% Define variables:
%  output:
%  cnt      -- Barycenter of the core. 
%
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
%   01/15/2018        0.9             hme
%

crv=CddCoreVertices(v);
scv=size(crv,1);
if scv==1
   cnt=crv;
else
   cnt=sum(crv)/scv;
end

function cnt=BaryCenter(clv)
% BaryCENTER computes the barycenter of the core.
%
% Usage: cnt=clv.BaryCenter()
% Define variables:
%  output:
%  cnt      -- Barycenter of the core. 
%
%  input:
%  clv      -- TuGame class object
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

crv=clv.CddCoreVertices();
scv=size(crv,1);
if scv==1
   cnt=crv;
else
   cnt=sum(crv)/scv;
end

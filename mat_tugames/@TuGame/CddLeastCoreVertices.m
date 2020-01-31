function LCV=CddLeastCoreVertices(clv,idx,tol)
% CDDLEASTCOREVERTICES computes the vertices of the least core of game v
% using CDDMEX.
%
%  Usage: LVC=clv.CddLeastCoreVertices(idx,tol)
%
% Define variables:
%  output:
%  Fields of LCV:
%
%  crv        -- Matrix of core vertices. Output is numeric.
%  crst       -- The core constraints.
%  vol        -- Least Core volume, if the core is full dimensional,
%                otherwise zero. 
%  P          -- Returns V- and H-representation (class Polyhedron)
%
%  epsv       -- Critical epsilon value of the least core.
%
%  input:
%  clv        -- TuGame class object.
%  idx        -- Specifies the projection plane onto R^3. Default is the empty set.
%                Will be computed internally.
%  tol        -- A positive tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/22/2014        0.5             hme
%

if nargin<3
 idx=[];   
 tol=10^6*eps;
else
  tol=10^6*eps;   
end

fmin=clv.CddLeastCore(tol);
v_eps=clv.streps_value(fmin);
[core_vert,crst,cr_vol,P]=CddCoreVertices(v_eps,idx,tol);
LCV=struct('crv',core_vert,'crst',crst,'vol',cr_vol,'P',P,'epsv',fmin);

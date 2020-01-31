function ALCV=CddAntiLeastCoreVertices(clv,idx,tol)
% CDDANTILEASTCOREVERTICES computes the vertices of the anti least core of game v
% using CDDMEX.
%
%  Usage: ALVC=clv.CddAntiLeastCoreVertices(tol)
%
% Define variables:
%  output:
%  Fields of ALCV:
%
%  crv        -- Matrix of anti core vertices. Output is numeric.
%  crst       -- The anti core constraints.
%  vol        -- Anti Least Core volume, if the core is full dimensional,
%                otherwise zero. 
%  P          -- Returns V- and H-representation (class Polyhedron)
%
%  epsv       -- Critical epsilon value of the anti least core.
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
%   08/09/2016        0.9             hme
%

if nargin<3
 idx=[];   
 tol=10^6*eps;
else
  tol=10^6*eps;   
end


fmin=clv.CddAntiLeastCore(tol);
v_eps=clv.streps_value(-fmin);
[core_vert,crst,cr_vol,P]=CddAntiCoreVertices(v_eps,idx,tol);
ALCV=struct('crv',core_vert,'crst',crst,'vol',cr_vol,'P',P,'epsv',fmin);

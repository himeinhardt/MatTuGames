function ALCV=CddAntiLeastCoreVertices(v,idx,tol)
% CDDANTILEASTCOREVERTICES computes the vertices of the anti least core of game v
% using CDDMEX.
%
%  Usage: ALVC=CddAntiLeastCoreVertices(v,tol)
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
%  v          -- A Tu-Game v of length 2^n-1. 
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

N=length(v);
[l1, n]=log2(N);

fmin=CddAntiLeastCore(v,tol);
v_eps=streps_value(v,-fmin);
[core_vert,crst,cr_vol,P]=CddAntiCoreVertices(v_eps,idx,tol);
ALCV=struct('crv',core_vert,'crst',crst,'vol',cr_vol,'P',P,'epsv',fmin);

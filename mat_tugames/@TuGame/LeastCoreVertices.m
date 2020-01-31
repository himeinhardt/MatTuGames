function [LCV,status]=LeastCoreVertices(clv,tol)
% LEASTCOREVERTICES computes the vertices of the least core of game v.
%
%  Usage: [crv,fmin,status]=clv.LeastCoreVertices(tol)
%
% Define variables:
%  output:
%  Fields of LCV:
%
%  crv        -- Matrix of core vertices. Output is numeric.
%  crst       -- The core constraints.
%  epsv       -- Critical epsilon value of the least core.
%  
%  Second Argument:
%  status     -- Returns 1 if the optimization terminated
%                successfully.
%
%  input:
%  clv        -- TuGame class object.
%  tol        -- A positive tolerance value. Its default value is set to 10^8*eps.
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
 tol=10^8*eps;
else
  tol=10^8*eps;   
end

[fmin,~,~,~,status]=clv.LeastCore(tol);
v_eps=clv.streps_value(fmin);
[core_vert,crst]=CoreVertices(v_eps);
LCV=struct('crv',core_vert,'crst',crst,'epsv',fmin);

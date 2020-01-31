function [LCV,status]=LeastCoreVertices(v,tol)
% LEASTCOREVERTICES computes the vertices of the least core of game v.
%
%  Usage: [crv,fmin,status]=LeastCoreVertices(v,tol)
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
%  v          -- A Tu-Game v of length 2^n-1. 
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

N=length(v);
[l1, n]=log2(N);

[fmin,~,~,~,status]=LeastCore(v,tol);
v_eps=streps_value(v,fmin);
[core_vert,crst]=CoreVertices(v_eps);
LCV=struct('crv',core_vert,'crst',crst,'epsv',fmin);

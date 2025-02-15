function ALsol=Anti_LorenzSol(v,tol)
% ANTI_LORENZSSOL determines the anti-Lorenz solution of game v. 
% Requires Matlab's Optimization Toolbox or cplexmex from 
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.5.1 and higher)
%
% Usage: ALsol=Anti_LorenzSol(v,tol)
%
% Define fields of the structure variable Lsol:
%
%  output:
%  ACp       -- Closest point of the anti-core to x, that is, the anti-Lorenz solution.
%  D         -- The distance between the points.
%  acr_vaild -- Indicates if ACp is in the anti-core (1 otherwise 0). 
%  resid     -- The residual.
%  ef        -- Exitflag.
%  lambda    -- Containing the Lagrange multipliers.
%  x         -- The reference point from which the distance to core should be drawn.
%              that is the equal distribution of the grand coalition.
%
%  input:
%  v            -- A Tu-Game v of length 2^n-1. 
%  tol          -- A positive tolerance value. Its default value is set to 10^6*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/20/2021        1.0             hme
%


if nargin<2
  tol=10^6*eps;
end


N=length(v);
[~, n]=log2(N);

if CddAntiCoreQ(v,tol)==0
   ALsol=inf(1,n);
   return;
end

x=ones(1,n)*v(N)/n;
ALsol=Anti_CPCore(v,x,tol);
AbQ=belongToAntiCoreQ(v,x,'rat',tol);
if AbQ
  ALsol.ACp=x;
  ALsol.AD=0;
  ALsol.resid=zeros(1,n);
end


function Lsol=LorenzSol(v,tol)
% LORENZSSOL determines the Lorenz solution of game v. 
% Requires Matlab's Optimization Toolbox or cplexmex from 
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.5.1 and higher)
%
% Usage: Lsol=LorenzSol(v,tol)
%
% Define fields of the structure variable Lsol:
%
%  output:
%  Cp       -- Closest point of the core to x, that is, the Lorenz solution.
%  D        -- The distance between the points.
%  cr_vaild -- Indicates if Cp is in the core (1 otherwise 0). 
%  resid    -- The residual.
%  ef       -- Exitflag.
%  lambda   -- Containing the Lagrange multipliers.
%  x        -- The reference point from which the distance to core should be drawn.
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
%   01/13/2019        1.0             hme
%


if nargin<2
  tol=10^6*eps;
end


N=length(v);
[~, n]=log2(N);

if CddCoreQ(v)==0
   Lsol=inf(1,n);
   return;
end

x=ones(1,n)*v(N)/n;
Lsol=CPCore(v,x,tol);
bQ=belongToCoreQ(v,x);
if bQ
  Lsol.Cp=x;
  Lsol.D=0;
  Lsol.resid=zeros(1,n);
end


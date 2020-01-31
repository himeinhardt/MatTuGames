function DR=DuttaRay(v,tol)
%DUTTARAY computes the Dutta-Ray solution for convex games.
% Requures Matlab's Optimization Toolbox or cplexmex from the URL
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.5.1 and higher) 
%
% Usage: DR=DuttaRay(v,tol)
% Define fields of the structure variable DR:
%  output:
%  Cp       -- Closest point of the core to x.
%  D        -- The distance between the points.
%  cr_vaild -- Indicates if Cp is in the core (1), otherwise false (0). 
%  resid    -- The residual.
%  ef       -- Exitflag.
%  lambda   -- Containing the Lagrange multipliers.
%  x        -- The reference point from which the distance to core should be drawn.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/14/2019        1.0             hme
%



if nargin<2
  tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
cvQ=convex_gameQ(v);
if cvQ==1
   x1=ones(1,n)*v(N)/n;
   bQ=belongToCoreQ(v,x1);
   DR=CPCore(v,x1,tol);
   if bQ
      DR.Cp=x1;
      DR.D=0;
      DR.resid=zeros(1,n);
   end
else
    DR=-inf(1,n);
end

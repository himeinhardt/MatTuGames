function LMQ=LorenzMaxCoreQ(v,x,tol)
% LorenzMaxCoreQ checks if x is Lorenz maximal in the core of game v, i.e., x is in the Lorenz set.
%
% Source: Hougaard, J.L., Peleg, B., Thorlund-Petersen, L., 2001. On the set of Lorenz-maximal
%         imputations in the core of a balanced game. Internat. J. Game Theory 30, 147–165.
%
% Usage: LMQ=LorenzMaxCoreQ(v,x,tol)
%
% Define variables:
%  output: field variables
%  LMQ      -- Returns 1 (true) whenever x is Lorenz maximal in the core, otherwise 0 (false).
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- An allocation x of length n.
%  tol      -- Tolerance value. Its default value is set to 10^10*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/20/2021        1.9.1           hme
%    
 
    
if nargin < 3
   tol=10^8*eps;
end    
try
   crQ=CddCoreQ(v);
catch
   crQ=coreQ(v);
end	

if crQ==1	
   bcQ=belongToCoreQ(v,x,'rat',tol);    
   crv=CddCoreVertices(v);
   s1=size(crv,1);
   lmQ=false(1,s1);
   dzQ=false(1,s1);
   for k=1:s1
       LD = LorenzDom(v,x,crv(k,:),tol);
       lmQ(k)=LD.ldQ;
       dmQ=abs(LD.dij.x)<=tol;
       dzQ(k)=all(dmQ(:));
   end
   LMQ=all(lmQ)  & all(dzQ) & bcQ;
else
   LMQ=false;	
end	

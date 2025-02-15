function ALMQ=Anti_LorenzMinACoreQ(clv,x,tol)
% Anti_LorenzMinCoreQ checks if x is Anti_Lorenz minimal in the anti-core of game v, i.e., x is in the Anti_Lorenz set.
%
% Source: Hougaard, J.L., Peleg, B., Thorlund-Petersen, L., 2001. On the set of Lorenz-maximal
%         imputations in the core of a balanced game. Internat. J. Game Theory 30, 147–165.
%
% Usage: ALMQ=clv.Anti_LorenzMinACoreQ(x,tol)
%
% Define variables:
%  output: field variables
%  ALMQ     -- Returns 1 (true) whenever x is anti-Lorenz minimal in the anti-core, otherwise 0 (false).
%
%  input:
%  clv      -- TuGame class object
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
%   11/24/2021        1.9.1           hme
%    
 
    
if nargin < 3
   tol=10^8*eps;
end    
try
   crQ=clv.CddAntiCoreQ(tol);
catch
   crQ=clv.anti_coreQ(tol);
end	

if crQ==1	
   bcQ=clv.belongToAntiCoreQ(x,'rat',tol);    
   crv=clv.CddAntiCoreVertices();
   s1=size(crv,1);
   lmQ=false(1,s1);
   dzQ=false(1,s1);
   for k=1:s1
       LD = clv.Anti_LorenzDom(x,crv(k,:),tol);
       lmQ(k)=LD.ldQ;
       dmQ=abs(LD.dij.x)<=tol;
       dzQ(k)=all(dmQ(:));
   end
   ALMQ=all(lmQ)  & all(dzQ) & bcQ;
else
   ALMQ=false;	
end	

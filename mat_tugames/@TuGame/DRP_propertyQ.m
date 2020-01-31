function DRP=DRP_propertyQ(clv,x,str,tol)
%DRP_PROPERTYQ checks wheter the solution x satisfies the dual replication property.
%
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          Sudhoelter (1997), The modified nucleolus: Properties and axiomatizations. International Journal of Game Theory, 26
%          (2):147–182, Jun 1997. ISSN 1432-1270. doi: 10.1007/BF01295846. URL https://doi.org/10.1007/BF01295846.
%
% Usage: DRP=clv.DRP_propertyQ(x,str,tol)
% Define variables:
%  output:
%  DRP     -- Returns true (1) whenever the solution fulfills dual replication property,
%              otherwise false (0).
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the dual extension of game v 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the dual extension of game v 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, the dual extension of game v 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the dual extension of game v.
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the dual extension of game v 
%               in accordance with modified pre-kernel solution.
%              'PMPRK' that is, the dual extension of game v 
%               in accordance with proper modified pre-kernel solution.
%              Default is 'MPRK'.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/01/2018        1.0             hme
%

if nargin<3
 tol=10^6*eps; % Change this value if the solution is not correct.
 str='MODIC';
elseif nargin<4
 tol=10^6*eps;
end

N=clv.tusize;
n=clv.tuplayers;

[~,v_ext]=clv.DualCover();
v_x=v_ext.d;
y=[x,x];

if strcmp(str,'SHAP')
     solQ=clv.ShapleyQ(x,tol);
     if solQ==1
        z=ShapleyValue(v_x);
        solQ_x=all(abs(y-z)<tol);
    else
      solQ_x=false;
    end
elseif strcmp(str,'PRN')
   solQ=clv.balancedCollectionQ(x,tol);
   if solQ==1
      solQ_x=balancedCollectionQ(v_x,y,tol);
   else
      solQ_x=false;
   end
elseif strcmp(str,'PRK')
   solQ=clv.PrekernelQ(x,tol);
   if solQ==1
      solQ_x=PrekernelQ(v_x,y,tol);
   else
      solQ_x=false;
   end
elseif strcmp(str,'MODIC')
   solQ=clv.modiclusQ(x,tol);
   if solQ==1
      solQ_x=modiclusQ(v_x,y,tol);
   else 
      solQ_x=false;
   end
elseif strcmp(str,'MPRK')
   solQ=clv.ModPrekernelQ(x,tol);
   if solQ==1
      solQ_x=ModPrekernelQ(v_x,y,tol);
   else
      solQ_x=false;
   end
elseif strcmp(str,'PMPRK')
   solQ=clv.PModPrekernelQ(x,tol);
   if solQ==1
      solQ_x=PModPrekernelQ(v_x,y,tol);
   else
      solQ_x=false;
   end
else 
   solQ=clv.ModPrekernelQ(x,tol);
   if solQ==1
      solQ_x=ModPrekernelQ(v_x,y,tol);
   else
      solQ_x=false;
   end
end 
DRP.propQ=solQ_x;
DRP.xQ=solQ;
DRP.y=y;
DRP.x=x;

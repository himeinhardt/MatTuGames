function SDFP=SDFP_propertyQ(clv,x,str,tol)
%SDFP_PROPERTYQ checks wheter the solution x satisfies a strong dual floor property.
%
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%  
%
% Usage: SDFP=clv.SDFP_propertyQ(x,str,tol)
%
% Define variables:
%
%  output:
%  DFP      -- Returns true (1) whenever the solution fulfills dual floor property,
%              otherwise false (0).
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'APRN' that is, the dual floor of game v 
%               in accordance with the anti-pre-nucleolus.
%              'APRK' that is, the dual floor of game v 
%               in accordance with anti-pre-kernel solution.
%              'SHAP' that is, the dual floor of game v 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the dual floor of game v.
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the dual floor of game v 
%               in accordance with modified pre-kernel solution.
%              'PMPRK' that is, the dual floor of game v 
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
%   07/16/2018        1.0             hme
%

N=clv.tusize;
n=clv.tuplayers;


if nargin<3
 tol=10^6*eps; % Change this value if the solution is not correct.
 str='MODIC';
elseif nargin<4
 tol=10^6*eps;
end


[~,ext_v]=clv.DualFloor();

y=[x,x];
v_x=ECFloorGame(ext_v.p,y);



if strcmp(str,'SHAP')
     solQ=clv.ShapleyQ(x,tol);
     if solQ==1
        z=ShapleyValue(v_x);
        y=[x,x];
        solQ_x=all(abs(y-z)<tol);
    else
      solQ_x=false;
    end
elseif strcmp(str,'APRN')
   solQ=clv.Anti_balancedCollectionQ(x,tol);
   if solQ==1
      solQ_x=Anti_balancedCollectionQ(v_x,y,tol);
   else
      solQ_x=false;
   end
elseif strcmp(str,'APRK')
   solQ=clv.Anti_PrekernelQ(x,tol);
   if solQ==1
      solQ_x=Anti_PrekernelQ(v_x,y,tol);
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
SDFP.propQ=solQ_x;
SDFP.xQ=solQ;
SDFP.y=y;
SDFP.x=x;

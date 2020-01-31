function DFP=DFP_propertyQ(clv,x,str,tol)
%DFP_PROPERTYQ checks wheter the solution x satisfies the dual floor property.
%
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%  
%
% Usage: DFP=clv.DFP_propertyQ(x,str,tol)
% Define variables:
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
%   05/06/2019        1.1             hme
%

if nargin<3
 tol=10^6*eps; % Change this value if the solution is not correct.
 str='MODIC';
elseif nargin<4
 tol=10^6*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
dc_v=clv.DualFloor();
if length(v)>1
   Mgc=clv.AllMarginalContributions();
   d=max(Mgc)-min(Mgc);
   t=6*sum(d);
else
   d=0;
   t=0;
end

v_x=streps_value(dc_v,t);
y=[x,x];

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
DFP.propQ=solQ_x;
DFP.xQ=solQ;
DFP.y=y;
DFP.x=x;

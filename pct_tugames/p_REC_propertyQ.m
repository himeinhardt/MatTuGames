function RECQ=p_REC_propertyQ(v,x,str,tol)
% P_REC_PROPERTYQ checks wheter the solution x satisfies reverse excess comparability using Matlab's PCT.
%
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: REC=p_REC_propertyQ(v,x,str,tol)
% Define variables:
%  output:
%  ECQ      -- Returns true (1) whenever the solution fulfills excess comparability,
%              otherwise false (0).
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'APRN' that is, the excess comparability floor 
%               in accordance with the anti pre-nucleolus.
%              'APRK' that is, the excess comparability floor
%               in accordance with anti pre-kernel solution.
%              'SHAP' that is, the excess comparability floor
%               in accordance with the Shapley Value.
%              'MODIC' that is, the excess comparability floor
%               equivalence in accordance with the modiclus.
%              'MAPRK' that is, the excess comparability floor game 
%               in accordance with modified anti-pre-kernel solution.
%              'PMAPRK' that is, the excess comparability floor game 
%               in accordance with proper modified anti-pre-kernel solution.
%              Default is 'MODIC'.
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

N=length(v);
[~, n]=log2(N);
v_x=p_ECFloorGame(v,x);

if strcmp(str,'SHAP')
     y=p_ShapleyValue(v_x);
elseif strcmp(str,'APRN')
   try
     y=cplex_Anti_PreNucl_llp(v_x);
   catch
     y=Anti_PreNucl_llp(v_x);
   end
elseif strcmp(str,'APRK')
     y=p_Anti_PreKernel(v_x,x);
elseif strcmp(str,'MODIC')
   dc_v=DualFloor(v_x);
   try 
     z=cplex_prenucl_llp(dc_v);
   catch
     z=PreNucl_llp(dc_v);
   end
   y=z(1:n);
elseif strcmp(str,'MAPRK')
     y=p_Anti_ModPreKernel(v_x,x);
elseif strcmp(str,'PMAPRK')
     y=p_Anti_PModPreKernel(v_x,x);
else
     y=p_Anti_ModPreKernel(v_x,x);
end 

df=all(abs(y-x)<tol);
RECQ.propQ=df;
RECQ.y=y;
RECQ.x=x;

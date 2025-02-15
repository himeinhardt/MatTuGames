function ECQ=p_EC_propertyQ(clv,x,str,tol)
% P_EC_PROPERTYQ checks whether the solution x satisfies excess comparability
% using Matlab's PCT.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          Sudhoelter (1997), The modified nucleolus: Properties and axiomatizations. International Journal of Game Theory, 26
%          (2):147–182, Jun 1997. ISSN 1432-1270. doi: 10.1007/BF01295846. URL https://doi.org/10.1007/BF01295846.
%
% Usage: EC=clv.p_EC_propertyQ(x,str,tol)
% Define variables:
%  output:
%  ECQ      -- Returns true (1) whenever the solution fulfills excess comparability,
%              otherwise false (0).
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods.
%              Permissible methods are:
%              'PRN' that is, the excess comparability cover game
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the excess comparability cover game
%               in accordance with pre-kernel solution.
%              'SHAP' that is, the excess comparability cover game
%               in accordance with the Shapley Value.
%              'MODIC' that is, the excess comparability cover game.
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the excess comparability cover game
%               in accordance with modified pre-kernel solution.
%              'PMPRK' that is, the excess comparability cover game
%               in accordance with proper modified pre-kernel solution.
%              Default is 'MPRK'.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional)
%
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
v_x=clv.p_ECCoverGame(x);

if strcmp(str,'SHAP')
     y=p_ShapleyValue(v_x);
elseif strcmp(str,'PRN')
   try
     y=msk_prenucl_llp(v_x); % MOSEK solver.
   catch
     y=PreNucl_llp(v_x);
   end
elseif strcmp(str,'PRK')
     y=p_PreKernel(v_x,x);
elseif strcmp(str,'MODIC')
   dc_v=DualCover(v_x);
   try 
     z=msk_prenucl_llp(dc_v); % MOSEK solver.
   catch
     z=PreNucl_llp(dc_v);
   end
   y=z(1:n);
elseif strcmp(str,'MPRK')
     y=p_ModPreKernel(v_x,x);
elseif strcmp(str,'PMPRK')
     y=p_PModPreKernel(v_x,x);
else
     y=p_ModPreKernel(v_x,x);
end 

df=all(abs(y-x)<tol);
ECQ.propQ=df;
ECQ.y=y;
ECQ.x=x;

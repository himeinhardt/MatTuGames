function ERGPQ=p_EC_RGP_Q(v,x,str,tol)
% P_ECRGPQ checks whether the solution x satisfies excess comparability for each reduced game
% using Matlab's PCT.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: ERGPQ=p_EC_RGP_Q(x,str,tol)
%
% Define variables:
%
%  output:
%  ECQ      -- Returns true (1) whenever each solution of the associated reduced game
%              fulfills excess comparability, otherwise false (0).
%  ecQ      -- Retruns a list of Boolean ones/zeros indicating for each indivudal
%              reduced game if EC is satisfied.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods.
%              Permissible methods are:
%              'PRN' that is, the maximum excess dual cover game
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the maximum excess dual cover game
%               in accordance with pre-kernel solution.
%              'SHAP' that is, the maximum excess dual cover game
%               in accordance with the Shapley Value.
%              'MODIC' that is, the maximum excess dual cover game.
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the maximum excess dual cover game
%               in accordance with modified pre-kernel solution.
%              'PMPRK' that is, the maximum excess dual cover game
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
%   06/16/2018        1.0             hme
%


if nargin<2
   N=length(v);
   tol=10^6*eps;
   str='MODIC';
   mnc_v=cplex_modiclus(v);
elseif nargin<3
   N=length(v);
   tol=10^6*eps;
   str='MODIC';
   mncQ=modiclusQ(v,x);
   if mncQ==1
      mnc_v=x;
   else
      warning('Sol:Wrn','Input vector is not the modiclus!');
      mnc_v=x;
      mnc_v=cplex_modiclus(v);
   end
elseif nargin < 4
   N=length(v);
   tol=10^6*eps;
   mncQ=modiclusQ(v,x);
   if mncQ==1
      mnc_v=x;
   else
      warning('Sol:Wrn','Input vector is not the modiclus!');
      mnc_v=x;
      %mnc_v=cplex_modiclus(v);
   end
else
   N=length(v);
   mncQ=modiclusQ(v,x);
   if mncQ==1
      mnc_v=x;
   else
      warning('Sol:Wrn','Input vector is not the modiclus!');
      mnc_v=x;
      %mnc_v=cplex_modiclus(v);
   end
end
[~, n]=log2(N);
vS=DM_Reduced_game(v,mnc_v);
ECQ=false;
ecQ=false(1,N-1);
k=1:n;
parfor S=1:N-1
    a=bitget(S,k)==1;
    try 
      x_S=cplex_modiclus(vS{1,S});
    catch
      x_S=Modiclus(vS{1,S});
    end  
    ECQ=EC_propertyQ(vS{1,S},x_S,str,tol);
    ecQ(S)=ECQ.propQ;
end

ECQ=all(ecQ);
ERGPQ.ECQ=ECQ;
ERGPQ.ecQ=ecQ;

function [SEDC SEDCGPQ]=p_Sedcons_propertyQ(v,x,str,tol)
% P_SEDCONS_propertyQ checks whether an imputation x satisfies the
% sedcons property (consistency) using Matlab's PCT.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: [SEDC SEDCGPQ]=p_Sedcons_propertyQ(v,x,str,tol)
%
% Define variables:
%
%  Output Fields
%  rgpQ     -- Returns 1 (true) whenever the SEDCONS is satisfied, 
%              otherwise 0 (false).
%  rgpq     -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  sedoncs  -- Returns a lits of reduced games for which the sedcons
%              property is satsified. An array of zeros and/or ones is returned.%  
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'APRN' that is, the Davis-Maschler anti-reduced game of the ECF game  
%               in accordance with the anti-pre-nucleolus.
%              'APRK' that is, the Davis-Maschler anti-reduced game of the ECF game 
%               in accordance with anti-pre-kernel solution.
%              'SHAP' that is, Hart-MasColell anti-reduced game of the ECF game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler anti-reduced game of the ECF game 
%               equivalence in accordance with the modiclus.
%              'MAPRK' that is, the Davis-Maschler anti-reduced game of the ECF game 
%               in accordance with the modified anti-pre-kernel solution.
%              'PMAPRK' that is, the Davis-Maschler anti-reduced game of the ECF game 
%               in accordance with the proper modified anti-pre-kernel solution.
%              'HMS_APK' that is, Hart-MasColell anti-reduced game of the ECF game 
%               in accordance with the anti-pre-kernel solution.
%              'HMS_APN' that is, Hart-MasColell anti-reduced game of the ECF game 
%               in accordance with the anti-pre-nucleous.
%              Default is 'MAPRK'
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 
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



if nargin<2
  x=p_Anti_ModPreKernel(v);
  n=length(x);
  tol=10^6*eps;
  str='MAPRK';
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='MAPRK';
elseif nargin<4
  n=length(x);
  tol=10^6*eps;
else
  n=length(x);
end

N=length(v);
S=1:N;
rgpq=false(1,N);
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2)==1;
impVec=cell(1,N);
rgpq_sol=cell(1,N);
sol=cell(1,N);
sedcons=false(1,N);

v_x=p_ECFloorGame(v,x);
%vS=cell(2,N);
if strcmp(str,'SHAP')
  vSa=p_HMS_AntiReduced_game(v_x,x,'SHAP');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_APK')
  vSa=p_HMS_AntiReduced_game(v_x,x,'APRK');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_APN')
  vSa=p_HMS_AntiReduced_game(v_x,x,'APRN');
  vS=vSa{:,1};
  clear vSa;  
else
  vSa=p_DM_AntiReduced_game(v_x,x);
  vS=vSa{:,1};
  clear vSa;
end


parfor k=1:N-1
 impVec{k}=x(PlyMat(k,:)); 
  if strcmp(str,'SHAP')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS.
   sol{k}=ShapleyValue(vS{k});
   rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
   SED=SED_propertyQ(vS{k},impVec{k});
   sedcons(k)=SED.propQ;
   rgpq(k)=all(rgpq_sol{k});
  elseif strcmp(str,'APRK')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS. To speed up computation, we use this code below for both, 
% the pre-nucleolus and and the pre-kernel. 
   SED=SED_propertyQ(vS{k},impVec{k});
   sedcons(k)=SED.propQ;
   rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'APRN')
   if length(vS{k})==1
     SED=SED_propertyQ(vS{k},impVec{k});
     sedcons(k)=SED.propQ;
     rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_AntiPreNucl(vS{k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{k}=Anti_PreNucl_llp(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     SED=SED_propertyQ(vS{k},impVec{k});
     sedcons(k)=SED.propQ;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'MODIC')
   if length(vS{k})==1
     SED=SED_propertyQ(vS{k},impVec{k});
     sedcons(k)=SED.propQ;
     rgpq(k)=modiclusQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_modiclus(vS{k}); % using cplex solver. 
     catch
       sol{k}=Modiclus(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     SED=SED_propertyQ(vS{k},impVec{k});
     sedcons(k)=SED.propQ;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'MAPRK')
   SED=SED_propertyQ(vS{k},impVec{k});
   sedcons(k)=SED.propQ;
   rgpq(k)=Anti_ModPrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'PMAPRK')
   SED=SED_propertyQ(vS{k},impVec{k});
   sedcons(k)=SED.propQ;
   rgpq(k)=Anti_PModPrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'HMS_APK')
   SED=SED_propertyQ(vS{k},impVec{k});
   sedcons(k)=SED.propQ;
   rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'HMS_APN')
   if length(vS{k})==1
     SED=SED_propertyQ(vS{k},impVec{k});
     sedcons(k)=SED.propQ;
     rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_AntiPreNucl(vS{k});
     catch
       sol{k}=Anti_PreNucl_llp(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     SED=SED_propertyQ(vS{k},impVec{k});
     sedcons(k)=SED.propQ;
     rgpq(k)=all(rgpq_sol{k});
   end
  else
     SED=SED_propertyQ(vS{k},impVec{k});
     sedcons(k)=SED.propQ;
     rgpq(k)=Anti_ModPrekernelQ(vS{k},impVec{k});   
  end
end

if strcmp(str,'SHAP')
   sol{N}=ShapleyValue(v_x);
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'APRK')
   sol{N}=Anti_PreKernel(v_x);
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;
   rgpq(N)=Anti_PrekernelQ(v_x,x);
elseif strcmp(str,'APRN')
   try
     sol{N}=cplex_AntiPreNucl_llp(v_x);
   catch
     sol{N}=Anti_PreNucl_llp(v_x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MODIC')
   try
     sol{N}=cplex_modiclus(v_x); % cplex_solver
   catch
     sol{N}=Modiclus(v_x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;   
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MAPRK')
   sol{N}=Anti_ModPreKernel(v_x);
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;
   rgpq(N)=Anti_ModPrekernelQ(v_x,x);
elseif strcmp(str,'PMAPRK')
   sol{N}=Anti_PModPreKernel(v_x);
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;
   rgpq(N)=Anti_PModPrekernelQ(v_x,x);
elseif strcmp(str,'HMS_APK')
  sol{N}=Anti_PreKernel(v_x);
  SED=SED_propertyQ(v_x,sol{N});
  sedcons(N)=SED.propQ;
  rgpq(N)=Anti_PrekernelQ(v_x,x);
elseif strcmp(str,'HMS_APN')
   try
     sol{N}=cplex_AntiPreNucl_llp(v_x);
   catch
     sol{N}=Anti_PreNucl_llp(v_x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;
   rgpq(N)=all(rgpq_sol{N});
else 
   sol{N}=Anti_ModPreKernel(v_x);
   SED=SED_propertyQ(v_x,sol{N});
   sedcons(N)=SED.propQ;
   rgpq(N)=Anti_ModPrekernelQ(v_x,x);
end
rgpQ=all(rgpq & sedcons);
%Formatting Output
if nargout>1
 SEDC=struct('sedconsQ',rgpQ,'rgpq',rgpq,'sedpropQ',sedcons);
 SEDCGPQ={'vS',vS,'impVec',impVec};
else
  SEDC=struct('sedconsQ',rgpQ,'rgpq',rgpq,'sedpropQ',sedcons);
end

function [LEDC LEDCGPQ]=Ledcons_propertyQ(v,x,str,tol)
% LEDCONS_propertyQ checks whether an imputation x satisfies the
% ledcons property (consistency).
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          Sudhoelter (1997), The modified nucleolus: Properties and axiomatizations. International Journal of Game Theory, 26
%          (2):147–182, Jun 1997. ISSN 1432-1270. doi: 10.1007/BF01295846. URL https://doi.org/10.1007/BF01295846.
%
% Usage: [LEDC LEDCGPQ]=Ledcons_propertyQ(v,x,str,tol)
%
% Define variables:
%
%  Output Fields
%  rgpQ     -- Returns 1 (true) whenever the LEDCONS is satisfied, 
%              otherwise 0 (false).
%  rgpq     -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and/or ones.
%  ledoncs  -- Returns a lits of reduced games for which the ledcons
%              property is satsified. An array of zeros and/or ones is returned.
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%
%  Input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the Davis-Maschler reduced game of the ECC game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Davis-Maschler reduced game of the ECC game
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game of the ECC game
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler reduced game of the ECC game
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the Davis-Maschler reduced game of the ECC game
%               in accordance with the modified pre-kernel solution.
%              'PMPRK' that is, the Davis-Maschler reduced game of the ECC game
%               in accordance with the proper modified pre-kernel solution.
%              'HMS_PK' that is, Hart-MasColell reduced game of the ECC game
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell reduced game of the ECC game
%               in accordance with the pre-nucleous.
%              Default is 'MPRK'.
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
  x=ModPreKernel(v);
  n=length(x);
  tol=10^6*eps;
  str='MPRK';
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='MPRK';
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
ledcons=false(1,N);

v_x=ECCoverGame(v,x);
%vS=cell(2,N);
if strcmp(str,'SHAP')
  vS=HMS_Reduced_game(v_x,x,'SHAP');
elseif strcmp(str,'HMS_PK')
  vS=HMS_Reduced_game(v_x,x,'PRK');
elseif strcmp(str,'HMS_PN')
  vS=HMS_Reduced_game(v_x,x,'PRN');
else
  vS=DM_Reduced_game(v_x,x);
end

for k=1:N-1
 impVec{1,k}=x(PlyMat(k,:)); 
  if strcmp(str,'SHAP')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS.
   sol{1,k}=ShapleyValue(vS{1,k});
   rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
   LED=LED_propertyQ(vS{1,k},impVec{1,k});
   ledcons(k)=LED.propQ;
   rgpq(k)=all(rgpq_sol{1,k});
  elseif strcmp(str,'PRK')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS. To speed up computation, we use this code below for both, 
% the pre-nucleolus and and the pre-kernel.
   LED=LED_propertyQ(vS{1,k},impVec{1,k});
   ledcons(k)=LED.propQ; 
   rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'PRN')
   if length(vS{1,k})==1
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=Prenucl(vS{1,k},impVec{1,k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{1,k}=PreNucl2(vS{1,k},impVec{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'MODIC')
   if length(vS{1,k})==1
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=modiclusQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=msk_modiclus(vS{1,k}); % using mosek modiclus function. 
     catch
       sol{1,k}=Modiclus(vS{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'MPRK')
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=ModPrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'PMPRK')
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=PModPrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'HMS_PK')
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'HMS_PN')
   if length(vS{1,k})==1
      LED=LED_propertyQ(vS{1,k},impVec{1,k});
      ledcons(k)=LED.propQ;
      rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=Prenucl(vS{1,k},impVec{1,k});
     catch
       sol{1,k}=PreNucl2(vS{1,k},impVec{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  else
     LED=LED_propertyQ(vS{1,k},impVec{1,k});
     ledcons(k)=LED.propQ;
     rgpq(k)=ModPrekernelQ(vS{1,k},impVec{1,k});   
  end
end

if strcmp(str,'SHAP')
   sol{N}=ShapleyValue(v_x);
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'PRK')
   sol{N}=PreKernel(v_x);
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=PrekernelQ(v_x,x);
elseif strcmp(str,'PRN')
   try
     sol{N}=Prenucl(v_x,x);
   catch
     sol{N}=PreNucl2(v_x,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MODIC')
   try
     sol{N}=msk_modiclus(v_x);
   catch
     sol{N}=Modiclus(v_x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;   
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MPRK')
   sol{N}=ModPreKernel(v_x);
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=ModPrekernelQ(v_x,x);
elseif strcmp(str,'PMPRK')
   sol{N}=PModPreKernel(v_x);
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=PModPrekernelQ(v_x,x);
elseif strcmp(str,'HMS_PK')
   sol{N}=PreKernel(v_x);
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=PrekernelQ(v_x,x);
elseif strcmp(str,'HMS_PN')
   try
     sol{N}=Prenucl(v_x,x);
   catch
     sol{N}=PreNucl2(v_x,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=all(rgpq_sol{N});
else 
   sol{N}=ModPreKernel(v_x);
   LED=LED_propertyQ(v_x,sol{N});
   ledcons(N)=LED.propQ;
   rgpq(N)=ModPrekernelQ(v_x,x);
end
rgpQ=all(rgpq & ledcons);
%Formatting Output
if nargout>1
 LEDC=struct('ledconsQ',rgpQ,'rgpq',rgpq,'ledpropQ',ledcons);
 LEDCGPQ={'vS',vS,'impVec',impVec};
else
  LEDC=struct('ledconsQ',rgpQ,'rgpq',rgpq,'ledpopQ',ledcons);
end

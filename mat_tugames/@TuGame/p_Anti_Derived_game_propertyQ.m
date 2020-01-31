function [ADGP ADGPC]=p_Anti_Derived_game_propertyQ(clv,x,str,tol)
% P_ANTI_DERIVED_GAME_PROPERTYQ checks whether an imputation x satisfies 
% a modified anti-derived game property (consistency) using Matlab's PCT.
%
% Usage: [ADGP ADGPC]=clv.Anti_Derived_game_propertyQ(x,str,tol)
% Define variables:
%  output: Fields
%  dgpQ     -- Returns 1 (true) whenever the ADGP is satisfied, 
%              otherwise 0 (false).
%  dgpq     -- Gives a precise list of anti-derived games for which the 
%              restriction of x on S is a solution of the anti-derived game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell anti-derived games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'APRN' that is, the Davis-Maschler anti-derived game 
%               in accordance with the anti-pre-nucleolus.
%              'APRK' that is, the Davis-Maschler anti-derived game 
%               in accordance with anti-pre-kernel solution.
%              'SHAP' that is, Hart-MasColell anti-derived game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler anti-derived game.
%               equivalence in accordance with the modiclus.
%              'PMAPRK' that is, the Davis-Maschler anti-derived game.
%               equivalence in accordance with the proper modified anti-pre-kernel.
%              'MAPRK' that is, the Davis-Maschler derived game.
%               equivalence in accordance with the modified anti-pre-kernel.
%              'HMS_APK' that is, Hart-MasColell anti-derived game 
%               in accordance with the anti-pre-kernel solution.
%              'HMS_APN' that is, Hart-MasColell anti-derived game 
%               in accordance with the anti-pre-nucleous.
%              'HMS_MODIC' that is, the Hart-MasColell anti-derived game
%               equivalence in accordance with the modiclus.
%              Default is 'MODIC'.
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
%   06/25/2018        1.0             hme
%                



if nargin<2
  x=clv.Modcilus()
  tol=10^6*eps;
  str='MODIC';
  n=clv.tuplayers;
elseif nargin<3
  tol=10^6*eps;
  str='MODIC';
  n=clv.tuplayers;
elseif nargin<4
  tol=10^6*eps;
  n=clv.tuplayers;
else
  n=clv.tuplayers;
end

N=clv.tusize;
S=1:N;
rgpq=false(1,N);
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2)==1;
impVec=cell(1,N);
rgpq_sol=cell(1,N);
sol=cell(1,N);

%vS=cell(2,N);
vS=cell(N-1,1);
if strcmp(str,'SHAP')
  vS=clv.p_HMS_Anti_Derived_game(x,'SHAP');
elseif strcmp(str,'HMS_APK')
  vS=clv.p_HMS_Anti_Derived_game(x,'APRK');
elseif strcmp(str,'HMS_APN')
  vS=clv.p_HMS_Anti_Derived_game(x,'APRN');
elseif strcmp(str,'HMS_MODIC')
  vS=clv.p_HMS_Anti_Derived_game(x,'MODIC');
else
  vS=clv.p_DM_Anti_Derived_game(x);
end


parfor k=1:N-1
 impVec{k}=x(PlyMat(k,:)); 
  if strcmp(str,'SHAP')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS.
   sol{k}=ShapleyValue(vS{k});
   rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
   rgpq(k)=all(rgpq_sol{k});
  elseif strcmp(str,'APRK')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS. To speed up computation, we use this code below for both, 
% the pre-nucleolus and and the pre-kernel. 
   rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'MAPRK')
   rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'PMAPRK')
   rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'APRN')
   if length(vS{k})==1
     rgpq(k)=Anti_balancedCollectionQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_AntiPreNucl_llp(vS{k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{k}=Anti_PreNucl_llp(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'MODIC')
   if length(vS{k})==1
     rgpq(k)=modiclusQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_AntiPreNucl_llp(vS{k}); % using cplex pre-nucleolus function.
     catch
       sol{k}=Anti_PreNucl_llp(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'HMS_APK')
   rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'HMS_APN')
   if length(vS{k})==1
     rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_AntiPreNucl_llp(vS{k});
     catch
       sol{k}=Anti_PreNucl_llp(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'HMS_MODIC')
   if length(vS{k})==1
     rgpq(k)=modiclusQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_AntiPreNucl_llp(vS{k});
     catch
       sol{k}=Anti_PreNucl_llp(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end   
  end
end

v=vS{1,N};

if strcmp(str,'SHAP')
   sol{N}=p_ShapleyValue(v);
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'APRK')
  rgpq(N)=Anti_PrekernelQ(v,x);
elseif strcmp(str,'MAPRK')
  rgpq(N)=PrekernelQ(v,x);
elseif strcmp(str,'PMAPRK')
  rgpq(N)=PrekernelQ(v,x);
elseif strcmp(str,'APRN')
   try
     sol{N}=cplex_AntiPreNucl_llp(v);
   catch
     sol{N}=Anti_PreNucl_llp(v); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MODIC')
   try
     sol{N}=cplex_AntiPreNucl_llp(v);
   catch
     sol{N}=Anti_PreNucl_llp(v); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_APK')
  rgpq(N)=Anti_PrekernelQ(v,x);
elseif strcmp(str,'HMS_APN')
   try
     sol{N}=cplex_AntiPreNucl_llp(v);
   catch
     sol{N}=Anti_PreNucl_llp(v,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_MODIC')
   try
     sol{N}=cplex_AntiPreNucl_llp(v);
   catch
     sol{N}=Anti_PreNucl_llp(v); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
end
rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 ADGP=struct('dgpQ',rgpQ,'dgpq',rgpq);
 ADGPC={'vS',vS,'impVec',impVec};
else
  ADGP=struct('dgpQ',rgpQ,'dgpq',rgpq);
end

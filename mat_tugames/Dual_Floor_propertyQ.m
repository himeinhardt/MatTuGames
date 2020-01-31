function [DGP DGPC]=Dual_Floor_propertyQ(v,x,str,tol)
% DERIVED_GAME_PROPERTYQ checks whether an imputation x satisfies 
% a modified anti-reduced game property (consistency).
%
% Usage: [DGP DGPC]=Dual_Floor_propertyQ(v,x,str,tol)
% Define variables:
%  output: Fields
%  dgpQ     -- Returns 1 (true) whenever the DGP is satisfied, 
%              otherwise 0 (false).
%  dgpq     -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'APRN' that is, the Davis-Maschler derived game 
%               in accordance with the anti-pre-nucleolus.
%              'APRK' that is, the Davis-Maschler anti-derived game 
%               in accordance with anti-pre-kernel solution.
%              'SHAP' that is, Hart-MasColell anti-derived game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler anti-derived game.
%               equivalence in accordance with the modiclus.
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
%   02/06/2018        0.9             hme
%                



if nargin<2
  x=Modcilus(v);
  n=length(x);
  tol=10^6*eps;
  str='MODIC';
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='MODIC';
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

%vS=cell(2,N);
if strcmp(str,'SHAP')
  vS=HMS_Anti_Derived_game(v,x,'SHAP');
elseif strcmp(str,'HMS_APK')
  vS=HMS_Anti_Derived_game(v,x,'APRK');
elseif strcmp(str,'HMS_APN')
  vS=HMS_Anti_Derived_game(v,x,'APRN');
elseif strcmp(str,'HMS_MODIC')
  vS=HMS_Anti_Derived_game(v,x,'MODIC');
else
  vS=Dual_Floor_game(v,x);
end


for k=1:N-1
 impVec{1,k}=x(PlyMat(k,:)); 
  if strcmp(str,'SHAP')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS.
   sol{1,k}=ShapleyValue(vS{1,k});
   rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
   rgpq(k)=all(rgpq_sol{1,k});
  elseif strcmp(str,'APRK')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS. To speed up computation, we use this code below for both, 
% the pre-nucleolus and and the pre-kernel. 
   rgpq(k)=Anti_PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'APRN')
   if length(vS{1,k})==1
     rgpq(k)=Anti_balancedCollectionQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=Anti_PreNucl(vS{1,k},impVec{1,k}); % using anti-pre-nucleolus function.
     catch
       sol{1,k}=cplex_AntiPreNucl_llp(vS{1,k},impVec{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'MODIC')
   if length(vS{1,k})==1
     rgpq(k)=modiclusQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=cplex_AntiPreNucl_llp(vS{1,k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{1,k}=Anti_PreNucl_llp(vS{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'HMS_APK')
   rgpq(k)=Anti_PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'HMS_APN')
   if length(vS{1,k})==1
     rgpq(k)=Anti_PrekernelQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=Anti_PreNucl(vS{1,k},impVec{1,k});
     catch
       sol{1,k}=cplex_AntiPreNucl(vS{1,k},impVec{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'HMS_MODIC')
   if length(vS{1,k})==1
     rgpq(k)=modiclusQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=cplex_AntiPreNucl_llp(vS{1,k});
     catch
       sol{1,k}=Anti_PreNucl_llp(vS{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end   
  end
end

if strcmp(str,'SHAP')
   sol{N}=ShapleyValue(v);
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'APRK')
  rgpq(N)=Anti_PrekernelQ(v,x);
elseif strcmp(str,'APRN')
   try
     sol{N}=Anti_PreNucl(v,x);
   catch
     sol{N}=cplex_AntiPreNucl(v,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MODIC')
   try
     sol{N}=cplex_modiclus(v);
   catch
     sol{N}=Modiclus(v); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_APK')
  rgpq(N)=Anti_PrekernelQ(v,x);
elseif strcmp(str,'HMS_APN')
   try
     sol{N}=Anti_PreNucl(v,x);
   catch
     sol{N}=cplex_AntiPreNucl(v,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_MODIC')
   try
     sol{N}=cplex_modiclus(v);
   catch
     sol{N}=Modiclus(v); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
end
rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 DGP=struct('dgpQ',rgpQ,'dgpq',rgpq);
 DGPC={'vS',vS,'impVec',impVec};
else
  DGP=struct('dgpQ',rgpQ,'dgpq',rgpq);
end

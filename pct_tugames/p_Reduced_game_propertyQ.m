function [RGP RGPC]=p_Reduced_game_propertyQ(v,x,str,tol)
% P_REDUCED_GAME_PROPERTYQ checks whether an imputation x satisfies the
% reduced game property (consistency).
%
% Usage: [RGP RGPC]=p_Reduced_game_propertyQ(v,x,str,tol)
% Define variables:
%  output: Fields
%  rgpQ     -- Returns 1 (true) whenever the RGP is satisfied, 
%              otherwise 0 (false).
%  rgpq     -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the Davis-Maschler reduced game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Davis-Maschler reduced game 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, Davis-Maschler reduced game
%               equivalence in accordance with the modiclus
%              'HMS_PK' that is, Hart-MasColell reduced game 
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell reduced game 
%               in accordance with the pre-nucleous.
%              'HMS_MODIC' that is, Hart/MasColell reduced game
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the Davis-Maschler reduced game 
%               equivalence in accordance with the modified pre-kernel.
%              'PMPRK' that is, the Davis-Maschler reduced game 
%               equivalence in accordance with the proper modified pre-kernel.
%              Default is 'PRK'.
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
%   05/21/2011        0.1 alpha        hme
%   06/29/2012        0.2 beta         hme
%   05/29/2013        0.3              hme
%   05/16/2014        0.5              hme
%   05/16/2014        0.5              hme
%   02/10/2018        0.9              hme
%   04/09/2018        1.0              hme
%                



if nargin<2
  x=p_PreKernel(v);
  n=length(x);
  tol=10^6*eps;
  str='PRK';
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='PRK';
elseif nargin<4
  n=length(x);
  tol=10^6*eps;
else
  n=length(x);
end

N=length(v);
S=1:N;
rgpq=zeros(1,N);
PlyMat=false(N,n);
parfor i = 1:n,
    PlyMat(:,i) = bitget(S,i)==1;
end

impVec=cell(1,N);
rgpq_sol=cell(1,N);
sol=cell(1,N);

vS=cell(N-1,1);
if strcmp(str,'SHAP')
  vSa=p_HMS_Reduced_game(v,x,'SHAP');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_PK')
  vSa=p_HMS_Reduced_game(v,x,'PRK');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_PN')
  vSa=p_HMS_Reduced_game(v,x,'PRN');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_MODIC')
  vSa=p_HMS_Reduced_game(v,x,'MODIC');
  vS=vSa{:,1};
  clear vSa;
else
  vSa=p_DM_Reduced_game(v,x);
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
   rgpq(k)=all(rgpq_sol{k});
  elseif strcmp(str,'PRK')
% Checks whether a solution x restricted to S is a solution of the
% reduced game vS. To speed up computation, we use this code below for both,
% the pre-nucleolus and and the pre-kernel.
   rgpq(k)=PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'MPRK')
   rgpq(k)=ModPrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'PMPRK')
   rgpq(k)=PModPrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'PRN')
   if length(vS{k})==1
     rgpq(k)=PrekernelQ(vS{k},impVec{k});
   else
     try
       sol{k}=Prenucl(vS{k},impVec{k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{k}=PreNucl2(vS{k},impVec{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'MODIC')
   if length(vS{k})==1
     rgpq(k)=modiclusQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_modiclus(vS{k}); % using cplex pre-nucleolus function. 
     catch
       sol{k}=Modiclus(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'HMS_PK')
   rgpq(k)=PrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'HMS_PN')
   if length(vS{k})==1
     rgpq(k)=PrekernelQ(vS{k},impVec{k});
   else
     try
       sol{k}=Prenucl(vS{k},impVec{k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{k}=PreNucl2(vS{k},impVec{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end
  elseif strcmp(str,'HMS_MODIC')
   if length(vS{k})==1
     rgpq(k)=modiclusQ(vS{k},impVec{k});
   else
     try  
       sol{k}=cplex_modiclus(vS{k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{k}=Modiclus(vS{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end
  else
     rgpq(k)=PrekernelQ(vS{k},impVec{k});
  end
end


if strcmp(str,'SHAP')
   sol{N}=p_ShapleyValue(v);
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'PRK')
  rgpq(N)=p_PrekernelQ(v,x);
elseif strcmp(str,'MPRK')
  rgpq(N)=p_ModPrekernelQ(v,x);
elseif strcmp(str,'PMPRK')
  rgpq(N)=p_PModPrekernelQ(v,x);
elseif strcmp(str,'PRN')
   try
     sol{N}=Prenucl(v,x); % using adjusted Derks pre-nucleolus function.
  catch
     sol{N}=PreNucl2(v,x); % use a third party solver instead!
   end 
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MODIC')
   try
     sol{N}=cplex_modiclus(v); % using adjusted Derks pre-nucleolus function.
   catch
     sol{N}=Modiclus(v); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_PK')
  rgpq(N)=p_PrekernelQ(v,x);
elseif strcmp(str,'HMS_PN')
   try
     sol{N}=Prenucl(v,x); % using adjusted Derks pre-nucleolus function.
   catch
     sol{N}=PreNucl2(v,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_MODIC')
   try
      sol{N}=cplex_modiclus(v); % using adjusted Derks pre-nucleolus function.
   catch
      sol{N}=Modiclus(v); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
else 
  rgpq(N)=p_PrekernelQ(v,x);
end


rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
 RGPC={'vS',vS,'impVec',impVec};
else
  RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
end

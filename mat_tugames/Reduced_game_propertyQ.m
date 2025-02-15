function [RGP RGPC]=Reduced_game_propertyQ(v,x,str,tol)
% REDUCED_GAME_PROPERTYQ checks whether an imputation x satisfies the
% reduced game property (consistency).
%
% Usage: [RGP RGPC]=Reduced_game_propertyQ(v,x,str,tol)
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
%              'POPK' that is, the Davis-Maschler reduced game 
%               in accordance with positive pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler reduced game 
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the Davis-Maschler reduced game 
%               equivalence in accordance with the modified pre-kernel.
%              'PMPRK' that is, the Davis-Maschler reduced game 
%               equivalence in accordance with the proper modified pre-kernel.
%              'HMS_MODIC' that is, Hart/MasColell reduced game
%               equivalence in accordance with the modiclus.
%              'HMS_PK' that is, Hart-MasColell reduced game 
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell reduced game 
%               in accordance with the pre-nucleous.
%              'HMS_CORE' that is, the Hart-MasColell reduced game 
%               in accordance with the core.
%              'CORE' that is, the Davis-Maschler reduced game 
%               in accordance with the core.
%              'LOR' that is, the Davis-Maschler reduced game 
%               in accordance with the Lorenz solution.
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
%   08/19/2010        0.1 beta        hme
%   06/17/2012        0.2 beta        hme
%   05/28/2013        0.3             hme
%   02/06/2018        0.9             hme
%   04/07/2018        1.0             hme
%   06/18/2020        1.9             hme
%   10/17/2021        1.9.1           hme
%                



if nargin<2
  x=PreKernel(v);
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
rgpq=false(1,N);
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2)==1;
impVec=cell(1,N);
rgpq_sol=cell(1,N);
sol=cell(1,N);

%vS=cell(2,N);
if strcmp(str,'SHAP')
  vS=HMS_Reduced_game(v,x,'SHAP');
elseif strcmp(str,'HMS_PK')
  vS=HMS_Reduced_game(v,x,'PRK');
elseif strcmp(str,'HMS_PN')
  vS=HMS_Reduced_game(v,x,'PRN');
elseif strcmp(str,'HMS_CORE')
  vS=HMS_Reduced_game(v,x,'CORE');
elseif strcmp(str,'HMS_MODIC')
  vS=HMS_Reduced_game(v,x,'MODIC');
else
  vS=DM_Reduced_game(v,x);
end


for k=1:N-1
 impVec{1,k}=x(PlyMat(k,:)); 
  if strcmp(str,'SHAP')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS.
   sol{1,k}=ShapleyValue(vS{1,k});
   rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
   rgpq(k)=all(rgpq_sol{1,k});
  elseif strcmp(str,'PRK')
% Checks whether a solution x restricted to S is a solution of the 
% reduced game vS. To speed up computation, we use this code below for both, 
% the pre-nucleolus and and the pre-kernel. 
   rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'POPK')
   rgpq(k)=positivePrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'MPRK')
   rgpq(k)=ModPrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'PMPRK')
   rgpq(k)=PModPrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'CORE')
    try
       crQ=CddCoreQ(vS{1,k});
    catch
       crQ=coreQ(vS{1,k});
    end  
    if crQ==1    
      rgpq(k)=belongToCoreQ(vS{1,k},impVec{1,k},'rat',tol);
    else
      rgpq(k)=false;
    end	    
  elseif strcmp(str,'HMS_CORE')
    try
       crQ=CddCoreQ(vS{1,k});
    catch
       crQ=coreQ(vS{1,k});
    end 
    if crQ==1    
       rgpq(k)=belongToCoreQ(vS{1,k},impVec{1,k},'rat',tol);
    else
       rgpq(k)=false;
    end	    
  elseif strcmp(str,'PRN')
   if length(vS{1,k})==1
     rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=Prenucl(vS{1,k},impVec{1,k}); % using adjusted Derks pre-nucleolus function.
     catch
       sol{1,k}=PreNucl2(vS{1,k},impVec{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'MODIC')
   if length(vS{1,k})==1
     rgpq(k)=modiclusQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=msk_modiclus(vS{1,k}); % using mosek solver to get the modiclus. 
     catch
       sol{1,k}=Modiclus(vS{1,k}); % using default solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'LOR')
   rgpq(k)=LorenzMaxCoreQ(vS{1,k},impVec{1,k},tol);
  elseif strcmp(str,'HMS_PK')
   rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'HMS_PN')
   if length(vS{1,k})==1
     rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
   else
     try
       sol{1,k}=Prenucl(vS{1,k},impVec{1,k});
     catch
       sol{1,k}=PreNucl2(vS{1,k},impVec{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end
  elseif strcmp(str,'HMS_MODIC')
    if length(vS{1,k})==1
       rgpq(k)=modiclusQ(vS{1,k},impVec{1,k});
    else
     try
       sol{k}=msk_modiclus(vS{1,k}); % using mosek solver to get the modiclus.
     catch
       sol{k}=Modiclus(vS{1,k}); % using default solver instead!
     end
   rgpq_sol{k}=abs(sol{1,k}-impVec{1,k})<tol;
   rgpq(k)=all(rgpq_sol{1,k});
   end     
  else
    rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k}); 
  end
end

if strcmp(str,'SHAP')
   sol{N}=ShapleyValue(v);
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'LOR')
   rgpq(N)=LorenzMaxCoreQ(v,x,tol);      
elseif strcmp(str,'PRK')
  rgpq(N)=PrekernelQ(v,x);
elseif strcmp(str,'POPK')
  rgpq(N)=positivePrekernelQ(v,x);
elseif strcmp(str,'MPRK')
  rgpq(N)=ModPrekernelQ(v,x);
elseif strcmp(str,'PMPRK')
  rgpq(N)=PModPrekernelQ(v,x);
elseif strcmp(str,'CORE')
    try
        crQ=CddCoreQ(v);
    catch
        crQ=coreQ(v);
    end  
    if crQ==1    
       rgpq(N)=belongToCoreQ(v,x,'rat',tol);
    else
       rgpq(N)=false;
    end	    
elseif strcmp(str,'HMS_CORE')
    try
       crQ=CddCoreQ(v);
    catch
       crQ=coreQ(v);
    end
    if crQ==1
       rgpq(N)=belongToCoreQ(v,x,'rat',tol);
    else
       rgpq(N)=false;
    end	  
elseif strcmp(str,'PRN')
   try
     sol{N}=Prenucl(v,x);
   catch
     sol{N}=PreNucl2(v,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MODIC')
   try
     sol{N}=msk_modiclus(v); % using mosek solver to get the modiclus.
   catch
     sol{N}=Modiclus(v); % using default solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_PK')
  rgpq(N)=PrekernelQ(v,x);
elseif strcmp(str,'HMS_PN')
   try
     sol{N}=Prenucl(v,x);
   catch
     sol{N}=PreNucl2(v,x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_MODIC')
   try
     sol{N}=msk_modiclus(v); % using mosek solver to get the modiclus.
   catch
     sol{N}=Modiclus(v); % using default solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
else 
  rgpq(N)=PrekernelQ(v,x);
end
rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
 RGPC={'vS',vS,'impVec',impVec};
else
  RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
end

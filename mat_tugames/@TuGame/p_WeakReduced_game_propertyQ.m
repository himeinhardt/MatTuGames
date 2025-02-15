function [WRGP WRGPC]=p_WeakReduced_game_propertyQ(clv,x,str,tol)
% P_WEAKREDUCED_GAME_PROPERTYQ checks whether an imputation x satisfies the
% weak reduced game property (consistency) using MATLAB's PCT.
%
% Source: %          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%
% Usage: [WRGP WRGPC]=clv.p_WeakReduced_game_propertyQ(x,str,tol)
%
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
%              'PRN' that is, the Davis-Maschler weak reduced game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Davis-Maschler weak reduced game 
%               in accordance with pre-kernel solution.
%              'POPK' that is, the Davis-Maschler weak reduced game 
%               in accordance with positive pre-kernel solution.
%              'SHAP' that is, Hart-MasColell weak reduced game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler weak reduced game 
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the Davis-Maschler weak reduced game 
%               equivalence in accordance with the modified pre-kernel.
%              'PMPRK' that is, the Davis-Maschler weak reduced game 
%               equivalence in accordance with the proper modified pre-kernel.
%              'HMS_PK' that is, Hart-MasColell weak reduced game 
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell weak reduced game 
%               in accordance with the pre-nucleous.
%              'HMS_CORE' that is, the Hart-MasColell reduced game 
%               in accordance with the core.
%              'CORE' that is, the Davis-Maschler reduced game 
%               in accordance with the core.
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
%   11/03/2021        1.9.1           hme
%                


N=clv.tusize;
n=clv.tuplayers;

if nargin<2
   if isa(clv,'TuSol')
      x=clv.tu_prk;
   elseif isa(clv,'p_TuSol')
      x=clv.tu_prk;
   else
      x=clv.PreKernel();
   end
   if isempty(x)
     x=clv.PreKernel();
   end
  tol=10^6*eps;
  str='PRK';
elseif nargin<3
  tol=10^6*eps;
  str='PRK';
elseif nargin<4
  tol=10^6*eps;
else
  tol=10^6*eps;
end

if strcmp(str,'SHAP')
  [vSa,sS,PlyMat]=clv.p_HMS_TwoReduced_game(x,'SHAP');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_PK')
  [vSa,sS,PlyMat]=clv.p_HMS_TwoReduced_game(x,'PRK');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_PN')
  [vSa,sS,PlyMat]=clv.p_HMS_TwoReduced_game(x,'PRN');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_CORE')
  [vSa,sS,PlyMat]=clv.p_HMS_TwoReduced_game(x,'CORE');
  vS=vSa{:,1};
  clear vSa;
else
  [vSa,sS,PlyMat]=clv.p_DM_TwoReduced_game(x);
  vS=vSa{:,1};
  clear vSa;
end

lS=length(sS);
nlS=lS+1;
rgpq=false(1,nlS);
impVec=cell(1,nlS);
rgpq_sol=cell(1,nlS);
sol=cell(1,nlS);

parfor k=1:lS
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
  elseif strcmp(str,'POPK')
   rgpq(k)=positivePrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'MPRK')
   rgpq(k)=ModPrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'PMPRK')
   rgpq(k)=PModPrekernelQ(vS{k},impVec{k});
  elseif strcmp(str,'CORE')
    try
       crQ=CddCoreQ(vS{k});
    catch
       crQ=coreQ(vS{k});
    end  
    if crQ==1    
      rgpq(k)=belongToCoreQ(vS{k},impVec{k},'rat',tol);
    else
      rgpq(k)=false;
    end	
  elseif strcmp(str,'HMS_CORE')
    try
      crQ=CddCoreQ(vS{k});
    catch
      crQ=coreQ(vS{k});
    end
    if crQ==1
      rgpq(k)=belongToCoreQ(vS{k},impVec{k},'rat',tol);
    else
     rgpq(k)=false;
    end	
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
       sol{k}=msk_modiclus(vS{k}); % using msk pre-nucleolus function. 
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
       sol{k}=Prenucl(vS{k},impVec{k});
     catch
       sol{k}=PreNucl2(vS{k},impVec{k}); % use a third party solver instead!
     end
     rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
     rgpq(k)=all(rgpq_sol{k});
   end  
  else
    rgpq(k)=PrekernelQ(vS{k},impVec{k}); 
  end
end

if strcmp(str,'SHAP')
   sol{nlS}=clv.ShapleyValue();
   rgpq_sol{nlS}=abs(sol{nlS}-x)<tol;
   rgpq(nlS)=all(rgpq_sol{nlS});
elseif strcmp(str,'PRK')
  rgpq(nlS)=clv.PrekernelQ(x);
elseif strcmp(str,'POPK')
  rgpq(nlS)=clv.positivePrekernelQ(x);
elseif strcmp(str,'MPRK')
  rgpq(nlS)=clv.ModPrekernelQ(x);
elseif strcmp(str,'PMPRK')
  rgpq(nlS)=clv.PModPrekernelQ(x);
elseif strcmp(str,'PRN')
   try
     sol{nlS}=clv.Prenucl(x);
   catch
     sol{nlS}=clv.PreNucl2(x); % use a third party solver instead!
   end
   rgpq_sol{nlS}=abs(sol{nlS}-x)<tol;
   rgpq(nlS)=all(rgpq_sol{nlS});
elseif strcmp(str,'MODIC')
   try
     sol{nlS}=clv.msk_modiclus();
   catch
     sol{nlS}=clv.Modiclus(); % use a third party solver instead!
   end
   rgpq_sol{nlS}=abs(sol{nlS}-x)<tol;
   rgpq(nlS)=all(rgpq_sol{nlS});
elseif strcmp(str,'HMS_PK')
  rgpq(nlS)=clv.PrekernelQ(x);
elseif strcmp(str,'HMS_PN')
   try
     sol{nlS}=clv.Prenucl(x);
   catch
     sol{nlS}=clv.PreNucl2(x); % use a third party solver instead!
   end
   rgpq_sol{nlS}=abs(sol{nlS}-x)<tol;
   rgpq(nlS)=all(rgpq_sol{nlS});
elseif strcmp(str,'CORE')
    try
       crQ=clv.CddCoreQ();
    catch
       crQ=clv.coreQ();
    end
    if crQ==1
       rgpq(nlS)=clv.belongToCoreQ(x,'rat',tol);
    else
      rgpq(nlS)=false;
    end
elseif strcmp(str,'HMS_CORE')
    try
       crQ=clv.CddCoreQ();
    catch
       crQ=clv.coreQ();
    end
    if crQ==1
       rgpq(nlS)=clv.belongToCoreQ(x,'rat',tol);
    else
      rgpq(nlS)=false;
    end
else 
  rgpq(nlS)=clv.PrekernelQ(x);
end
rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 WRGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
 WRGPC={'vS',vS,'impVec',impVec};
else
  WRGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
end

function [DGP DGPC]=p_Derived_game_propertyQ(clv,x,str,tol)
% P_DERIVED_GAME_PROPERTYQ checks whether an imputation x satisfies a
% modifed reduced game property (consistency).
%
% Source: Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: [DGP DGPC]=clv.p_Derived_game_propertyQ(x,str,tol)
% Define variables:
%  output: Fields
%  dgpQ     -- Returns 1 (true) whenever the RGP is satisfied, 
%              otherwise 0 (false).
%  dgpq     -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the Davis-Maschler derived game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Davis-Maschler derived game 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell derived game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler derived game.
%               equivalence in accordance with the modiclus.
%              'PMPRK' that is, the Davis-Maschler derived game.
%               equivalence in accordance with the proper modified pre-kernel.
%              'MPRK' that is, the Davis-Maschler derived game.
%               equivalence in accordance with the modified pre-kernel.
%              'HMS_PK' that is, Hart-MasColell derived game 
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell derived game 
%               in accordance with the pre-nucleous.
%              'HMS_MODIC' that is, the Hart-MasColell derived game
%               equivalence in accordance with the modiclus.
%              Default is 'MODIC'.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional 
%              

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/01/2018        1.0              hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

if nargin<2
   x=clv.Modiclus();
  tol=10^6*eps;
  str='MODIC';
elseif nargin<3
  tol=10^6*eps;
  str='MODIC';
elseif nargin<4
  tol=10^6*eps;
else
  tol=10^6*eps;
end


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
  vS=clv.p_HMS_Derived_game(x,'SHAP');
elseif strcmp(str,'HMS_PK')
  vS=clv.p_HMS_Derived_game(x,'PRK');
elseif strcmp(str,'HMS_PN')
  vS=clv.p_HMS_Derived_game(x,'PRN');
elseif strcmp(str,'HMS_MODIC')
  vS=clv.p_HMS_Derived_game(x,'MODIC');
else
  vS=clv.p_DM_Derived_game(x);
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
   rgpq(k)=PrekernelQ(vS{1,k},impVec{k});
  elseif strcmp(str,'PMPRK')
   rgpq(k)=PrekernelQ(vS{1,k},impVec{k});
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
   if length(vS{1,k})==1
     rgpq(k)=modiclusQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_prenucl_llp(vS{k}); % using cplex pre-nucleolus function. 
     catch
       sol{k}=PreNucl_llp(vS{k}); % use a third party solver instead!
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
       sol{k}=cplex_prenucl_llp(vS{k});
     catch
       sol{k}=PreNucl_llp(vS{k}); % use a third party solver instead!
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
elseif strcmp(str,'PRK')
  rgpq(N)=p_PrekernelQ(v,x);
elseif strcmp(str,'MPRK')
  rgpq(N)=PrekernelQ(v,x);
elseif strcmp(str,'PMPRK')
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
     sol{N}=cplex_prenucl(v);
   catch
     sol{N}=PreNucl(v); % use a third party solver instead!
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
     sol{N}=cplex_prenucl(v);
   catch
     sol{N}=PreNucl(v); % use a third party solver instead!
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

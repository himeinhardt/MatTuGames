function [DGP DGPC]=Derived_game_propertyQ(clv,x,str,tol)
% DERIVED_GAME_PROPERTYQ checks whether an imputation x satisfies 
% a modified reduced game property (consistency).
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: [DGP DGPC]=clv.Derived_game_propertyQ(x,str,tol)
% Define variables:
%  output: Fields
%  dgpQ     -- Returns 1 (true) whenever the DGP is satisfied, 
%              otherwise 0 (false).
%  dgpq     -- Gives a precise list of derived games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell derived games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%
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
rgpq=false(1,N);
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2)==1;
impVec=cell(1,N);
rgpq_sol=cell(1,N);
sol=cell(1,N);

%vS=cell(2,N);
if strcmp(str,'SHAP')
  vS=clv.HMS_Derived_game(x,'SHAP');
elseif strcmp(str,'HMS_PK')
  vS=clv.HMS_Derived_game(x,'PRK');
elseif strcmp(str,'HMS_PN')
  vS=clv.HMS_Derived_game(x,'PRN');
elseif strcmp(str,'HMS_MODIC')
  vS=clv.HMS_Derived_game(x,'MODIC');
else
  vS=clv.DM_Derived_game(x);
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
  elseif strcmp(str,'MPRK')
   rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'PMPRK')
   rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k});
  elseif strcmp(str,'PRN')
   if length(vS{1,k})==1
     rgpq(k)=balancedCollectionQ(vS{1,k},impVec{1,k});
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
       sol{1,k}=cplex_prenucl_llp(vS{1,k}); % using cplex pre-nucleolus function. 
     catch
       sol{1,k}=PreNucl_llp(vS{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end
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
       sol{1,k}=cplex_prenucl_llp(vS{1,k});
     catch
       sol{1,k}=PreNucl_llp(vS{1,k}); % use a third party solver instead!
     end
     rgpq_sol{1,k}=abs(sol{1,k}-impVec{1,k})<tol;
     rgpq(k)=all(rgpq_sol{1,k});
   end   
  end
end

v=vS{1,N};

if strcmp(str,'SHAP')
   sol{N}=clv.ShapleyValue();
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'PRK')
  rgpq(N)=clv.PrekernelQ(x);
elseif strcmp(str,'MPRK')
  rgpq(N)=PrekernelQ(v,x);
elseif strcmp(str,'PMPRK')
  rgpq(N)=PrekernelQ(v,x);
elseif strcmp(str,'PRN')
   try
     sol{N}=clv.Prenucl(x);
   catch
     sol{N}=clv.PreNucl2(x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'MODIC')
   try
     sol{N}=clv.cplex_modiclus();
   catch
     sol{N}=clv.Modiclus(); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_PK')
  rgpq(N)=clv.PrekernelQ(x);
elseif strcmp(str,'HMS_PN')
   try
     sol{N}=clv.Prenucl(x);
   catch
     sol{N}=clv.PreNucl2(x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_MODIC')
   try
     sol{N}=clv.cplex_prenucl();
   catch
     sol{N}=clv.PreNucl(); % use a third party solver instead!
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

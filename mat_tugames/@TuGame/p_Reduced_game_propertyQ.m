function [RGP RGPC]=p_Reduced_game_propertyQ(clv,x,str,tol)
% P_REDUCED_GAME_PROPERTYQ checks whether an imputation x satisfies the
% reduced game property (consistency).
%
% Usage: [RGP RGPC]=clv.p_Reduced_game_propertyQ(x,str,tol)
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
%
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the Davis-Maschler reduced game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Davis-Maschler reduced game 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game 
%               in accordance with the Shapley Value.
%              'HMS_PK' that is, Hart-MasColell reduced game 
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell reduced game 
%               in accordance with the pre-nucleous.
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
%   05/30/2013        0.3              hme
%   05/16/2014        0.5              hme
%                

N=clv.tusize;
n=clv.tuplayers;

if nargin<2
   if isa(clv,'TuSol')
      x=clv.tu_prk;
   elseif isa(clv,'p_TuSol')
      x=clv.tu_prk;
   else
      x=clv.p_PreKernel();
   end
   if isempty(x)
     x=clv.p_PreKernel();
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
  vSa=clv.p_HMS_Reduced_game(x,'SHAP');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_PK')
  vSa=clv.p_HMS_Reduced_game(x,'PRK');
  vS=vSa{:,1};
  clear vSa;
elseif strcmp(str,'HMS_PN')
  vSa=clv.p_HMS_Reduced_game(x,'PRN');
  vS=vSa{:,1};
  clear vSa;
else
  vSa=clv.p_DM_Reduced_game(x);
  vS=vSa{:,1};
  clear vSa;
end

parfor k=1:N-1
 impVec{k}=x(PlyMat(k,:));
  if strcmp(str,'SHAP')
% Checks whether a solution x restricted to S is a solution of the
% reduced game vS.
   sol{k}=ShapleyValue(vS{1,k});
   rgpq_sol{k}=abs(sol{k}-impVec{k})<tol;
   rgpq(k)=all(rgpq_sol{k});
  elseif strcmp(str,'PRK')
% Checks whether a solution x restricted to S is a solution of the
% reduced game vS. To speed up computation, we use this code below for both,
% the pre-nucleolus and and the pre-kernel.
   rgpq(k)=PrekernelQ(vS{k},impVec{k});
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
  end
end


if strcmp(str,'SHAP')
   sol{N}=p_ShapleyValue(clv);
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'PRK')
  rgpq(N)=clv.p_PrekernelQ(x);
elseif strcmp(str,'PRN')
   try
     sol{N}=clv.Prenucl(x); % using adjusted Derks pre-nucleolus function.
   catch
     sol{N}=clv.PreNucl2(x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'HMS_PK')
  rgpq(N)=clv.p_PrekernelQ(x);
elseif strcmp(str,'HMS_PN')
   try
      sol{N}=clv.Prenucl(x); % using adjusted Derks pre-nucleolus function.
   catch
      sol{N}=clv.PreNucl2(x); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
end


rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
 RGPC={'vS',vS,'impVec',impVec};
else
  RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
end

function [RGP RGPC]=MaxConsistencyQ(v,x,str,tol)
% MAXCONSISTENCYQ checks whether an imputation x satisfies maximal
% consistency.
%
% Usage: [RGP RGPC]=MaxConsistencyQ(v,x,str,tol)
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
%              'HMS_PK' that is, Hart-MasColell reduced game 
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell reduced game 
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
%   06/10/2020        1.9             hme
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
N1=N-1;
rgpq=false(1,N1);
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2)==1;
impVec=cell(1,N1);
rgpq_sol=cell(1,N1);
sol=cell(1,N1);

%vS=cell(2,N);
if strcmp(str,'SHAP')
  vS=HMS_Reduced_game(v,x,'SHAP');
elseif strcmp(str,'HMS_PK')
  vS=HMS_Reduced_game(v,x,'PRK');
elseif strcmp(str,'HMS_PN')
  vS=HMS_Reduced_game(v,x,'PRN');
elseif strcmp(str,'HMS_CORE')
  vS=HMS_Reduced_game(v,x,'CORE');
else
  vS=DM_Reduced_game(v,x);
end


for k=1:N1
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
   rgpq(k)=belongToCoreQ(vS{1,k},impVec{1,k});  
  elseif strcmp(str,'HMS_CORE')
   rgpq(k)=belongToCoreQ(vS{1,k},impVec{1,k});	  
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
       sol{1,k}=cplex_modiclus(vS{1,k}); % using cplex pre-nucleolus function. 
     catch
       sol{1,k}=Modiclus(vS{1,k}); % use a third party solver instead!
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
  else
    rgpq(k)=PrekernelQ(vS{1,k},impVec{1,k}); 
  end
end


rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
 RGPC={'vS',vS,'impVec',impVec};
else
  RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
end

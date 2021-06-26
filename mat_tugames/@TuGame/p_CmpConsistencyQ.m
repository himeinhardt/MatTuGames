function [RGP RGPC]=p_CmpConsistencyQ(clv,x,str,tol)
% P_CMPCONSISTENCYQ checks whether an imputation x satisfies the complement
% consistency using Matlab's PCT.
%
% Source: Moulin H (1985) The separability axiom and equal sharing methods. J of Econ Theory 36:120-148
%
%
% Usage: [RGP RGPC]=clv.p_CmpConsistencyQ(x,str,tol)
% Define variables:
%  output: Fields
%  rgpQ     -- Returns 1 (true) whenever the complement consistency is satisfied, 
%              otherwise 0 (false).
%  rgpq     -- Gives a precise list of complement reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All complement reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the complement reduced game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the complement reduced game 
%               in accordance with pre-kernel solution.
%              'POPK' that is, the complement reduced game 
%               in accordance with positive pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the complement reduced game 
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the complement reduced game 
%               equivalence in accordance with the modified pre-kernel.
%              'PMPRK' that is, the complement reduced game 
%               equivalence in accordance with the proper modified pre-kernel.
%              'HMS_PK' that is, Hart-MasColell reduced game 
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell reduced game 
%               in accordance with the pre-nucleous.
%              'CORE' that is, the complement reduced game 
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
%   06/13/2020        1.9             hme
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

S=1:N;
N1=N-1;
rgpq=false(1,N1);
PlyMat=false(N,n);
parfor i = 1:n,
    PlyMat(:,i) = bitget(S,i)==1;
end
impVec=cell(1,N1);
rgpq_sol=cell(1,N1);
sol=cell(1,N1);

%vS=cell(2,N);
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
  vSa=clv.p_Complement_Reduced_game(x);
  vS=vSa{:,1};
  clear vSa;
end


for k=1:N1
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
   rgpq(k)=belongToCoreQ(vS{k},impVec{k});   
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


rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
 RGPC={'vS',vS,'impVec',impVec};
else
  RGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
end

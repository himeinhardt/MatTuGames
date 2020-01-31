function [ARGP ARGPC]=p_AntiReduced_game_propertyQ(clv,x,str,tol)
% P_ANTIREDUCED_GAME_PROPERTYQ checks whether an imputation x satisfies the
% anti-reduced game property (consistency) using Matlab's PCT.
%
% Usage: [ARGP ARGPC]=clv.p_AntiReduced_game_propertyQ(x,str,tol)
% Define variables:
%  Output Fields
%  rgpQ     -- Returns 1 (true) whenever the ARGP is satisfied, 
%              otherwise 0 (false).
%  rgpq     -- Gives a precise list of anti-reduced games for which the 
%              restriction of x on S is a solution of the anti-reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell anti-reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of x on all S.
%
%  Input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'APRN' that is, the Davis-Maschler anti-reduced game 
%               in accordance with the anti-pre-nucleolus.
%              'APRK' that is, the Davis-Maschler anti-reduced game 
%               in accordance with anti-pre-kernel solution.
%              'SHAP' that is, Hart-MasColell anti-reduced game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, checking covariance with strategic
%               equivalence in accordance with the modiclus
%              'HMS_APK' that is, Hart-MasColell anti-reduced game 
%               in accordance with the anti-pre-kernel solution.
%              'HMS_APN' that is, Hart-MasColell anti-reduced game 
%               in accordance with the anti-pre-nucleous.
%              Default is 'APRK'.
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



if nargin<2
  x=clv.Anti_PreKernel();
  tol=10^6*eps;
  str='APRK';
  n=clv.tuplayers;
elseif nargin<3
  tol=10^6*eps;
  str='APRK';
  n=clv.tuplayers;
elseif nargin<4
  tol=10^6*eps;
  n=clv.tuplayers;
else
  n=clv.tuplayers;
end

v=clv.tuvalues;
N=clv.tusize;

S=1:N;
rgpq=false(1,N);
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2)==1;
impVec=cell(1,N);
rgpq_sol=cell(1,N);
sol=cell(1,N);

%vS=cell(2,N);
if strcmp(str,'SHAP')
  vSa=clv.p_HMS_AntiReduced_game(x,'SHAP');
elseif strcmp(str,'HMS_APK')
  vSa=clv.p_HMS_AntiReduced_game(x,'APRK');
elseif strcmp(str,'HMS_APN')
  vSa=clv.p_HMS_AntiReduced_game(x,'APRN');
else
  vSa=clv.p_DM_AntiReduced_game(x);
end
vS=vSa{:,1};
clear vSa;


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
  elseif strcmp(str,'APRN')
   if length(vS{k})==1
     rgpq(k)=Anti_PrekernelQ(vS{k},impVec{k});
   else
     try
       sol{k}=cplex_AntiPreNucl_llp(vS{k}); % cplex solver.
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
       sol{k}=cplex_modiclus(vS{k}); % cplex solver 
     catch
       sol{k}=Modiclus(vS{k}); % use a third party solver instead!
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
  end
end

if strcmp(str,'SHAP')
   sol{N}=clv.p_ShapleyValue();
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
elseif strcmp(str,'APRK')
  rgpq(N)=clv.Anti_PrekernelQ(x);
elseif strcmp(str,'APRN')
   try
     sol{N}=clv.cplex_AntiPreNucl_llp();
   catch
     sol{N}=clv.Anti_PreNucl_llp(); % use a third party solver instead!
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
elseif strcmp(str,'HMS_APK')
  rgpq(N)=clv.Anti_PrekernelQ(x);
elseif strcmp(str,'HMS_APN')
   try
     sol{N}=clv.cplex_AntiPreNucl_llp();
   catch
     sol{N}=clv.Anti_PreNucl_llp(); % use a third party solver instead!
   end
   rgpq_sol{N}=abs(sol{N}-x)<tol;
   rgpq(N)=all(rgpq_sol{N});
end
rgpQ=all(rgpq);
%Formatting Output
if nargout>1
 ARGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
 ARGPC={'vS',vS,'impVec',impVec};
else
  ARGP=struct('rgpQ',rgpQ,'rgpq',rgpq);
end

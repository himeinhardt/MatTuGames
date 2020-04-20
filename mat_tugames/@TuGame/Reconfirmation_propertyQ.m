function [RCP RCPC]=Reconfirmation_propertyQ(clv,x,str,tol)
% RECONFIRMATION_PROPERTYQ checks whether an imputation x satisfies the RCP.
%
% Usage: [RCP RCPC]=clv.Reconfirmation_propertyQ(x,str,tol)
%
% Define variables:
%  output: Fields
%  RCPQ     -- Returns 1 (true) whenever the RCP is satisfied, 
%              otherwise 0 (false).
%  rcpq     -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  vS_sol   -- Returns a vector of solutions x_s for all reduced games vS.
%  vS_y     -- Returns a vector of extended solutions x_s to x_N for 
%              all reduced games vS.
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
%               in accordance with the pre-nucleolus.
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
%   05/29/2013        0.3             hme
%   04/01/2020        1.9             hme
%                

v=clv.tuvalues;
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
rcpQ=false(1,N);
PlyMat=false(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k)==1;end 
rcpq=cell(1,N);
vS_sol=cell(1,N);
vS_y=cell(1,N);

if strcmp(str,'SHAP')
  vS=clv.HMS_Reduced_game(x,'SHAP');
elseif strcmp(str,'HMS_PK')
  vS=clv.HMS_Reduced_game(x,'PRK');
elseif strcmp(str,'HMS_PN')
  vS=clv.HMS_Reduced_game(x,'PRN');
else
  vS=clv.DM_Reduced_game(x);
end

rS=cell(N-1);
for k=1:N-1
  vS_y{1,k}=x;
  if strcmp(str,'PRK')
    vS_sol{1,k}=PreKernel(vS{1,k}); % solution y restricted to S.
    rS{k}=PlyMat(k,:);
    vS_y{1,k}(rS{k})=vS_sol{1,k}; % extension to (y,x_N\S).
    rcpq{k}=clv.PrekernelQ(vS_y{1,k});
  elseif strcmp(str,'PRN')
    if length(vS{1,k})==1
      vS_sol{1,k}=PreKernel(vS{1,k}); % solution y restricted to S.
      rS{k}=PlyMat(k,:);
      vS_y{1,k}(rS{k})=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    else
      try
        vS_sol{1,k}=cplex_prenucl_mod4(vS{1,k}); % solution y restricted to S.
      catch
         vS_sol{1,k}=PreNucl(vS{1,k}); % use a third party solver instead!
      end
      rS{k}=PlyMat(k,:);
      vS_y{1,k}(rS{k})=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    end
  elseif strcmp(str,'SHAP')
    vS_sol{1,k}=ShapleyValue(vS{1,k}); % solution y restricted to S.
    rS{k}=PlyMat(k,:);
    vS_y{1,k}(rS{k})=vS_sol{1,k};     % extension to (y,x_N\S).
    rcpq{k}=abs(vS_y{1,k}-x)<tol;
  elseif strcmp(str,'HMS_PN')
    if length(vS{1,k})==1
      vS_sol{1,k}=PreKernel(vS{1,k}); % solution y restricted to S.
      rS{k}=PlyMat(k,:);
      vS_y{1,k}(rS{k})=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    else
      try
        vS_sol{1,k}=cplex_prenucl_mod4(vS{1,k}); % solution y restricted to S.
      catch
        vS_sol{1,k}=PreNucl(vS{1,k}); % use a third party solver instead!
      end
      rS{k}=PlyMat(k,:);
      vS_y{1,k}(rS{k})=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    end
 else % default
    vS_sol{1,k}=PreKernel(vS{1,k}); % solution y restricted to S.
    rS{k}=PlyMat(k,:);
    vS_y{1,k}(rS{k})=vS_sol{1,k}; % extension to (y,x_N\S)
    rcpq{k}=clv.PrekernelQ(vS_y{1,k});
  end
 rcpQ(k)=all(rcpq{k});
end

if strcmp(str,'PRK')
  vS_sol{1,N}=clv.PreKernel(x);
  rcpq{N}=clv.PrekernelQ(x);
  vS_y{1,N}=vS_sol{1,N};
elseif strcmp(str,'PRN')
  try
     vS_sol{1,N}=clv.cplex_prenucl_mod4();
  catch
     vS_sol{1,N}=clv.PreNucl();
  end
  rcpq{N}=abs(vS_sol{1,N}-x)<tol;
  vS_y{1,N}=vS_sol{1,N};
elseif strcmp(str,'SHAP')
  vS_sol{1,N}=clv.ShapleyValue();
  rcpq{N}=abs(vS_sol{1,N}-x)<tol;
  vS_y{1,N}=vS_sol{1,N};
else
  vS_sol{1,N}=clv.PreKernel(x);
  rcpq{N}=clv.PrekernelQ(x);
  vS_y{1,N}=vS_sol{1,N};
end
rcpQ(N)=all(rcpq{N});
RCPQ=all(rcpQ);
%Formatting Output
if nargout>1
 RCP=struct('RCPQ',RCPQ,'rcpQ',rcpQ);
 RCPC={'vS',vS,'vS_sol',vS_sol,'vS_y',vS_y};
else
  RCP=struct('RCPQ',RCPQ,'rcpQ',rcpQ);
end

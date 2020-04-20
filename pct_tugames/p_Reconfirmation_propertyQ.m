function [RCP RCPC]=p_Reconfirmation_propertyQ(v,x,str,tol)
% P_RECONFIRMATION_PROPERTYQ checks whether an imputation x satisfies the RCP.
%
% Usage: [RCP RCPC]=p_Reconfirmation_propertyQ(v,x,str,tol)
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
%   05/21/2011        0.1 alpha        hme
%   06/29/2012        0.2 beta         hme
%   05/27/2013        0.3              hme
%   05/16/2014        0.5              hme
%   04/01/2020        1.9             hme
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
rcpQ=zeros(1,N);
PlyMat=false(N,n);
parfor i = 1:n, PlyMat(:,i) = bitget(S,i)==1; end

rcpq=cell(1,N);
vS_sol=cell(1,N);
vS_y=cell(1,N);

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
else
  vSa=p_DM_Reduced_game(v,x);
  vS=vSa{:,1};
  clear vSa;
end

parfor k=1:N-1
  vS_y{1,k}=x;
  if strcmp(str,'PRK')
    vS_sol{1,k}=PreKernel(vS{k}); % solution y restricted to S.
    rSk=PlyMat(k,:);
    vS_y{1,k}(rSk)=vS_sol{1,k}; % extension to (y,x_N\S).
    rcpq{k}=PrekernelQ(v,vS_y{1,k});
  elseif strcmp(str,'PRN')
    if length(vS{k})==1
      vS_sol{1,k}=PreKernel(vS{k}); % solution y restricted to S.
      rSk=PlyMat(k,:);
      vS_y{1,k}(rSk)=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    else
      try
        vS_sol{1,k}=cplex_prenucl_mod4(vS{k}); % solution y restricted to S.
      catch
        vS_sol{1,k}=PreNucl(vS{k}); % use a third party solver instead!
      end
      rSk=PlyMat(k,:);
      vS_y{1,k}(rSk)=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    end
  elseif strcmp(str,'SHAP')
    vS_sol{1,k}=ShapleyValue(vS{k}); % solution y restricted to S.
    rSk=PlyMat(k,:);
    vS_y{1,k}(rSk)=vS_sol{1,k};     % extension to (y,x_N\S).
    rcpq{k}=abs(vS_y{1,k}-x)<tol;
  elseif strcmp(str,'HMS_PN')
    if length(vS{k})==1
      vS_sol{1,k}=PreKernel(vS{k}); % solution y restricted to S.
      rSk=PlyMat(k,:);
      vS_y{1,k}(rSk)=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    else
      try  
        vS_sol{1,k}=cplex_prenucl_mod4(vS{k}); % solution y restricted to S.
      catch
        vS_sol{1,k}=PreNucl(vS{k}); % use a third party solver instead!
      end
      rSk=PlyMat(k,:);
      vS_y{1,k}(rSk)=vS_sol{1,k}; % extension to (y,x_N\S).
      rcpq{k}=abs(vS_y{1,k}-x)<tol;
    end
 else % default
    vS_sol{1,k}=PreKernel(vS{k}); % solution y restricted to S.
    rSk=PlyMat(k,:);
    vS_y{1,k}(rSk)=vS_sol{1,k}; % extension to (y,x_N\S)
    rcpq{k}=PrekernelQ(v,vS_y{1,k});
  end
 rcpQ(k)=all(rcpq{k});
end

if strcmp(str,'PRK')
  vS_sol{1,N}=p_PreKernel(v,x);
  rcpq{N}=p_PrekernelQ(v,x);
  vS_y{1,N}=vS_sol{1,N};
elseif strcmp(str,'PRN')
  try 
    vS_sol{1,N}=cplex_prenucl_mod4(v);
  catch
    vS_sol{1,N}=PreNucl(v); % use a third party solver instead!
  end
  rcpq{N}=abs(vS_sol{1,N}-x)<tol;
  vS_y{1,N}=vS_sol{1,N};
elseif strcmp(str,'SHAP')
  vS_sol{1,N}=ShapleyValue(v);
  rcpq{N}=abs(vS_sol{1,N}-x)<tol;
  vS_y{1,N}=vS_sol{1,N};
else
  vS_sol{1,N}=p_PreKernel(v,x);
  crgpq{N}=p_PrekernelQ(v,x);
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

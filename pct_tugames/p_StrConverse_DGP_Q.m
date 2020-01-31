function [SCDGP SCDGPC]=p_StrConverse_DGP_Q(v,x,str,tol)
% P_STRCONVERSE_DGP_Q checks whether an imputation x satisfies the strong converse derived game 
% property (SCDPG) using Matlab's PCT, that is, the converse DGP (converse modified consistency).
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%
% Usage [SCDGP SCDGPC]=p_StrConverse_DGP_Q(v,x,str,tol)
%
% Define variables:
%
%  output: Fields
%  Q        -- Returns 1 (true) whenever the strong CDGP is satisfied, 
%              otherwise 0 (false).
%  scdgpQ   -- Gives a precise list of derived games for which the 
%              restriction of x on S is a solution of the derived game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell derived games on S 
%              with two players at x.
%  sV_x     -- Returns a vector of extended solutions x_s to x_N for 
%              all derived games vS.              
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'MPRK' that is, the Davis-Maschler derived game 
%               in accordance with the modified pre-kernel solution.
%              'PMPRK' that is, the Davis-Maschler derived game 
%               in accordance with the proper modified pre-kernel solution.
%              'SHAP' that is, Hart-MasColell derived game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler derived game 
%               equivalence in accordance with the modiclus
%              'PRN' that is, the Davis-Maschler derived game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Davis-Maschler derived game 
%               in accordance with pre-kernel solution.
%              Default is 'MPRK'.
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
%   03/25/2018        1.0             hme
%                

if nargin<2
  x=p_ModPreKernel(v);
  n=length(x);
  tol=10^6*eps;
  str='MPRK';
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='MPRK';
elseif nargin<4
  n=length(x);
  tol=10^6*eps;
else
  n=length(x);
end

N=length(v);
S=1:N;
crgpQ=false(1,N);
PlyMat=false(N,n);
parfor k=1:n, PlyMat(:,k) = bitget(S,k)==1;end

crgpq=cell(1,N);
impVec=cell(1,N);
sV_sol=cell(1,N);
sV_x=cell(1,N);
rS=cell(1,N);

vS=cell(2,N);
if strcmp(str,'SHAP')
  vS=p_HMS_Derived_game(v,x,'SHAP');
elseif strcmp(str,'HMS_MPK')
  vS=p_HMS_Derived_game(v,x,'PRK');
elseif strcmp(str,'HMS_PMPK')
  vS=p_HMS_Derived_game(v,x,'PRK');
else
  vS=p_DM_Derived_game(v,x);
end


parfor k=1:N-1
 sV_x{k}=x;
 impVec{k}=x(logical(PlyMat(k,:)));
  if strcmp(str,'MPRK')
    sV_sol{k}=PreKernel(vS{k},impVec{k});  % solution x restricted to S.
    rS{k}=PlyMat(k,:);
    sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
  elseif strcmp(str,'PMPRK')
    sV_sol{k}=PreKernel(vS{k},impVec{k});  % solution x restricted to S.
    rS{k}=PlyMat(k,:);
    sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
  elseif strcmp(str,'MODIC')
    if length(impVec{k})==1
       sV_sol{k}=Modiclus(vS{k},impVec{k});
       rS{k}=PlyMat(k,:);
       sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
    else
       try
         sV_sol{k}=cplex_prenucl(vS{k}); % solution x restricted to S.
       catch
         sV_sol{k}=PreNucl(vS{k}); % use a third party solver instead!
       end
       rS{k}=PlyMat(k,:);
       sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
    end
  elseif strcmp(str,'PRK')
    sV_sol{k}=PreKernel(vS{k},impVec{k});  % solution x restricted to S.
    rS{k}=PlyMat(k,:);
    sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
  elseif strcmp(str,'PRN')
    if length(impVec{k})==1
       sV_sol{k}=PreKernel(vS{k},impVec{k});
       rS{k}=PlyMat(k,:);
       sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
    else
       try
         sV_sol{k}=cplex_prenucl_llp(vS{k}); % solution x restricted to S.
       catch
         sV_sol{k}=PreNucl_llp(vS{k}); % use a third party solver instead!
       end
       rS{k}=PlyMat(k,:);
       sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
    end
  elseif strcmp(str,'SHAP')
    sV_sol{k}=ShapleyValue(vS{k}); % solution x restricted to S.
    rS{k}=PlyMat(k,:);
    sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
  else
    sV_sol{k}=PreKernel(vS{k},impVec{k});  % solution x restricted to S.
    rS{k}=PlyMat(k,:);
    sV_x{k}(rS{k})=sV_sol{k}; % extension to (x,x_N\S).
  end
  crgpq{k}=abs(sV_x{k}-x)<tol;
  crgpQ(k)=all(crgpq{k});
end

if strcmp(str,'MPRK')
  sV_sol{N}=p_ModPreKernel(v,x);
elseif strcmp(str,'PMPRK')
  sV_sol{N}=p_PModPreKernel(v,x);
elseif strcmp(str,'MODIC')
  try
    sV_sol{N}=cplex_modiclus(v);
  catch
    sV_sol{N}=Modilcus(v); % use a third party solver instead!
  end
elseif strcmp(str,'PRK')
  sV_sol{N}=p_PreKernel(v,x);
elseif strcmp(str,'PRN')
  try
    sV_sol{N}=cplex_prenucl_llp(v);
  catch
    sV_sol{N}=PreNucl_llp(v); % use a third party solver instead!
  end
elseif strcmp(str,'SHAP')
  sV_sol{N}=p_ShapleyValue(v);
else
  sV_sol{N}=p_ModPreKernel(v,x);
end
crgpq{N}=abs(sV_sol{1,N}-x)<tol;
crgpQ(N)=all(crgpq{N});
CrgpQ=all(crgpQ);
%Formatting Output
if nargout>1
 SCDGP=struct('Q',CrgpQ,'scdgpQ',crgpQ);
 SCDGPC={'vS',vS,'sV_x',sV_x};
else
  SCDGP=struct('Q',CrgpQ,'scdgpQQ',crgpQ);
end

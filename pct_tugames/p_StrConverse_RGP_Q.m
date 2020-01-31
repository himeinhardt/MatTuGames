function [CRGP CRGPC]=p_StrConverse_RGP_Q(v,x,str,tol)
% P_StrConverse_RGP_Q checks whether an imputation x satisfies the strong CRGP, 
% that is, the strong converse reduced game property (converse consistency). 
% Note that, strong CRGP is stronger than the usual definition of CRGP, 
% since it checks whether CRGP holds for all S subsets of N.
% Using Matlab's PCT.
%
% Usage: [CRGP CRGPC]=p_StrConverse_RGP_Q(v,x,str,tol)
% Define variables:
%  output: Fields
%  CrgpQ    -- Returns 1 (true) whenever the strong CRGP is satisfied, 
%              otherwise 0 (false).
%  crgpQ    -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  impVec   -- Returns a vector of restrictions of solution x on S for all S. 
%  sV_sol   -- Returns a vector of solutions x_s for all reduced games vS.
%  sV_x     -- Returns a vector of extended solutions x_s to x_N for 
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
%  tol=10^6*eps;
  n=length(x);
end

N=length(v);
S=1:N;
crgpQ=zeros(1,N);
PlyMat=false(N,n);
parfor i = 1:n, PlyMat(:,i) = bitget(S,i)==1; end

crgpq=cell(1,N);
impVec=cell(1,N);
sV_sol=cell(1,N);
sV_x=cell(1,N);

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
 sV_x{1,k}=x;
 impVec{1,k}=x(PlyMat(k,:));
  if strcmp(str,'PRK')
    sV_sol{1,k}=PreKernel(vS{k},impVec{1,k});  % solution x restricted to S.
    rSk=PlyMat(k,:);
    sV_x{1,k}(rSk)=sV_sol{1,k}; % extension to (x,x_N\S).
  elseif strcmp(str,'PRN')
    if length(impVec{1,k})==1
       sV_sol{1,k}=PreKernel(vS{k},impVec{1,k});
       rSk=PlyMat(k,:);
       sV_x{1,k}(rSk)=sV_sol{1,k}; % extension to (x,x_N\S).
    else
       try
         sV_sol{1,k}=Prenucl(vS{k}); % solution x restricted to S.
       catch
         sV_sol{1,k}=PreNucl(vS{k}); % use a third party solver instead!
       end
       rSk=PlyMat(k,:);
       sV_x{1,k}(rSk)=sV_sol{1,k}; % extension to (x,x_N\S).
    end
  elseif strcmp(str,'SHAP')
    sV_sol{1,k}=ShapleyValue(vS{k}); % solution x restricted to S.
    rSk=PlyMat(k,:);
    sV_x{1,k}(rSk)=sV_sol{1,k}; % extension to (x,x_N\S).
  else
    sV_sol{1,k}=PreKernel(vS{k},impVec{1,k});  % solution x restricted to S.
    rSk=PlyMat(k,:);
    sV_x{1,k}(rSk)=sV_sol{1,k}; % extension to (x,x_N\S).
  end
  crgpq{k}=abs(sV_x{1,k}-x)<tol;
  crgpQ(k)=all(crgpq{k});
end

if strcmp(str,'PRK')
  sV_sol{1,N}=p_PreKernel(v,x);
elseif strcmp(str,'PRN')
  try
    sV_sol{1,N}=Prenucl(v);
  catch
    sV_sol{1,N}=PreNucl(v); % use a third party solver instead!
  end
elseif strcmp(str,'SHAP')
  sV_sol{1,N}=p_ShapleyValue(v);
else
  sV_sol{1,N}=p_PreKernel(v,x);
end
crgpq{N}=abs(sV_sol{1,N}-x)<tol;
crgpQ(N)=all(crgpq{N});
CrgpQ=all(crgpQ);
%Formatting Output
if nargout>1
 CRGP=struct('CrgpQ',CrgpQ,'crgpQ',crgpQ);
 CRGPC={'vS',vS,'sV_x',sV_x};
else
  CRGP=struct('CrgpQ',CrgpQ,'crgpQ',crgpQ);
end

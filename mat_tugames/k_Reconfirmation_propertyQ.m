function [kRCP kRCPC]=k_Reconfirmation_propertyQ(v,x,K,str,tol)
% k_RECONFIRMATION_PROPERTYQ checks whether an imputation x satisfies the RCP.
%
%
% Usage: [kRCP kRCPC]=k_Reconfirmation_propertyQ(v,x,K,str,tol)
% Define variables:
%  output: Fields
%  kRCPQ     -- Returns 1 (true) whenever the k-RCP is satisfied, 
%              otherwise 0 (false).
%  krcpQ     -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS        -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%              of size 2<=|S|<=K.
%  vS_sol    -- Returns a vector of solutions x_s for all reduced games vS.
%  vS_y      -- Returns a vector of extended solutions x_s to x_N for 
%              all reduced games vS.
%  input:
%  v         -- A Tu-Game v of length 2^n-1. 
%  x         -- payoff vector of size(1,n). Must be efficient.
%  str       -- A string that defines different Methods. 
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
%  tol       -- Tolerance value. By default, it is set to 10^6*eps.
%               (optional) 
%              


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/01/2010        0.1 beta        hme
%   06/19/2010        0.2 beta        hme
%   05/27/2013        0.3             hme
%                



if nargin<2
  x=PreKernel(v);
  n=length(x);
  tol=10^6*eps;
  str='PRK';
  K=2;
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='PRK';
  K=2;
elseif nargin<4
  n=length(x);
  tol=10^6*eps;
  str='PRK';
elseif nargin<5
  n=length(x);
  tol=10^6*eps;
else
  n=length(x);
end

if K<2
 error('K must be an integer greater than 2!');
elseif K>n
 error('K must be an integer equal to or smaller than n!');
 else
end


N=length(v);
S=1:N;
PlyMat=false(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k)==1;end

sumPM=PlyMat*ones(n,1);
kRCPQ=0;

vS=cell(K);

for int=1:K;
  slcl2=sumPM==int;
  clK{int}=S(slcl2);
  PlyMatK{int}=PlyMat(clK{int},:);
  lgtK(int)=length(clK{int});

 for k=1:lgtK(int)
  if strcmp(str,'SHAP')
    vS{int,k}=HMS_RedGame(v,x,clK{int}(k));
  elseif strcmp(str,'HMS_PK')
    vS{int,k}=HMS_RedGame(v,x,clK{int}(k));
  elseif strcmp(str,'HMS_PN')
    vS{int,k}=HMS_RedGame(v,x,clK{int}(k));
  else
    vS{int,k}=RedGame(v,x,clK{int}(k));
  end
 end
end

clm=max(lgtK);
vS_sol=cell(K,clm);
vS_y=cell(K,clm);
rcpq=cell(K,clm);


for int=1:K
 for k=1:lgtK(int)
  vS_y{int,k}=x;
   if strcmp(str,'PRK')
     vS_sol{int,k}=PreKernel(vS{int,k}); % solution y restricted to S.
     rS{int,k}=PlyMatK{int}(k,:);
     vS_y{int,k}(rS{int,k})=vS_sol{int,k}; % extension to (y,x_N\S).
     rcpq{int,k}=PrekernelQ(v,vS_y{int,k});
   elseif strcmp(str,'PRN')
      if length(vS{int,k})==1
        vS_sol{int,k}=PreKernel(vS{int,k}); % solution y restricted to S.
        rS{int,k}=PlyMatK{int}(k,:);
        vS_y{int,k}(rS{int,k})=vS_sol{int,k}; % extension to (y,x_N\S).
        rcpq{int,k}=abs(vS_y{int,k}-x)<tol;
      else
        try
          vS_sol{int,k}=Prenucl(vS{int,k}); % solution y restricted to S.
        catch
          vS_sol{int,k}=PreNucl(vS{int,k}); % use a third party solver instead!
        end
        rS{int,k}=PlyMatK{int}(k,:);
        vS_y{int,k}(rS{int,k})=vS_sol{int,k}; % extension to (y,x_N\S).
        rcpq{int,k}=abs(vS_y{int,k}-x)<tol;
      end
   elseif strcmp(str,'SHAP')
     vS_sol{int,k}=ShapleyValue(vS{int,k}); % solution y restricted to S.
     rS{int,k}=PlyMatK{int}(k,:);
     vS_y{int,k}(rS{int,k})=vS_sol{int,k};     % extension to (y,x_N\S).
     rcpq{int,k}=abs(vS_y{int,k}-x)<tol;
   elseif strcmp(str,'HMS_PN')
     if length(vS{int,k})==1
        vS_sol{int,k}=PreKernel(vS{int,k}); % solution y restricted to S.
        rS{int,k}=PlyMatK{int}(k,:);
        vS_y{int,k}(rS{int,k})=vS_sol{int,k}; % extension to (y,x_N\S).
        rcpq{int,k}=abs(vS_y{int,k}-x)<tol;
      else
        try
           vS_sol{int,k}=Prenucl(vS{int,k}); % solution y restricted to S.
        catch
           vS_sol{int,k}=PreNucl(vS{int,k}); % use a third party solver instead!
        end
        rS{int,k}=PlyMatK{int}(k,:);
        vS_y{int,k}(rS{int,k})=vS_sol{int,k}; % extension to (y,x_N\S).
        rcpq{int,k}=abs(vS_y{int,k}-x)<tol;
      end
   else % default
     vS_sol{int,k}=PreKernel(vS{int,k}); % solution y restricted to S.
     rS{int,k}=PlyMatK{int}(k,:);
     vS_y{int,k}(rS{int,k})=vS_sol{int,k}; % extension to (y,x_N\S)
     rcpq{int,k}=PrekernelQ(v,vS_y{int,k});
   end
   rcpQ{1,int}(k)=all(rcpq{int,k});
 end
krcpQ(int)=all(rcpQ{1,int});
end

kRCPQ=all(all(krcpQ));
%Formatting Output
if nargout>1
 kRCP=struct('kRCPQ',kRCPQ,'krcpQ',krcpQ);
 kRCPC={'vS',vS,'vS_sol',vS_sol,'vS_y',vS_y};
else
  kRCP=struct('kRCPQ',kRCPQ,'krcpQ',krcpQ);    
end

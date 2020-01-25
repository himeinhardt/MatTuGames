function [kRCP kRCPC]=p_k_Reconfirmation_propertyQ(clv,x,K,str,tol)
% k_RECONFIRMATION_PROPERTYQ checks whether an imputation x satisfies the RCP.
% Using Matlab's PCT.
%
%
% Usage: [kRCP kRCPC]=p_k_Reconfirmation_propertyQ(clv,x,K,str,tol)
%
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
%
%  input:
%  clv       -- TuGame class object.
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
%   05/30/2013        0.3             hme
%                

N=clv.tusize;
n=clv.tuplayers;

if nargin<2
   if isa(clv,'TuSol')
      x=clv.tu_prk;
   elseif isa(clv,'p_TuSol')
      x=clv.tu_prk;
   else
      x=p_PreKernel(clv);
   end
   if isempty(x)
     x=p_PreKernel(clv);
   end
  tol=10^6*eps;
  str='PRK';
  K=2;
elseif nargin<3
  tol=10^6*eps;
  str='PRK';
  K=2;
elseif nargin<4
  tol=10^6*eps;
  str='PRK';
elseif nargin<5
  tol=10^6*eps;
else
  tol=10^6*eps;
end

if K<2
 error('K must be an integer greater than 2!');
elseif K>n
 error('K must be an integer equal to or smaller than n!');
 else
end


S=1:N;
PlyMat=false(N,n);
parfor i = 1:n, PlyMat(:,i) = bitget(S,i)==1; end


sumPM=PlyMat*ones(n,1);
kRCPQ=0;

vS=cell(K,1);


parfor int=1:K;
  slcl2=sumPM==int;
  clK{int}=S(slcl2);
  PlyMatK{int}=PlyMat(clK{int},:);
  lgtK(int)=length(clK{int});

  vSK=cell(1,lgtK(int));

 for k=1:lgtK(int)
  if strcmp(str,'SHAP')
    vSK{k}=HMS_RedGame(clv,x,clK{int}(k));
  elseif strcmp(str,'HMS_PK')
    vSK{k}=HMS_RedGame(clv,x,clK{int}(k));
  elseif strcmp(str,'HMS_PN')
    vSK{k}=HMS_RedGame(clv,x,clK{int}(k));
  else
    vSK{k}=RedGame(clv,x,clK{int}(k));
  end
 end
vS{int}=vSK;
end

clm=max(lgtK);
rcpQ=cell(1,K);
KvS_sol=cell(K,1);
KvS_y=cell(K,1);

parfor int=1:K
  rcpq=cell(1,lgtK(int));
  vS_sol=cell(1,lgtK(int));
  vS_y=cell(1,lgtK(int));
 for k=1:lgtK(int)
  vS_y{k}=x;
   if strcmp(str,'PRK')
     vk=cell2mat(vS{int}(k));
     vS_sol{k}=PreKernel(vk); % solution y restricted to S.
     rSk=PlyMatK{int}(k,:);
     vS_y{k}(rSk)=vS_sol{k}; % extension to (y,x_N\S).
     rcpq{k}=PrekernelQ(clv,vS_y{k});
   elseif strcmp(str,'PRN')
      vk=cell2mat(vS{int}(k));     
      if length(vk)==1
        vS_sol{k}=PreKernel(vk); % solution y restricted to S.
        rSk=PlyMatK{int}(k,:);
        vS_y{k}(rSk)=vS_sol{k}; % extension to (y,x_N\S).
        rcpq{k}=abs(vS_y{k}-x)<tol;
      else
        vk=cell2mat(vS{int}(k));
        try
          vS_sol{k}=Prenucl(vk); % solution y restricted to S.
        catch
          vS_sol{k}=PreNucl(vk); % use a third party solver instead!
        end
        rSk=PlyMatK{int}(k,:);
        vS_y{k}(rSk)=vS_sol{k}; % extension to (y,x_N\S).
        rcpq{k}=abs(vS_y{k}-x)<tol;
      end
   elseif strcmp(str,'SHAP')
     vk=cell2mat(vS{int}(k));
     vS_sol{k}=ShapleyValue(vk); % solution y restricted to S.
     rSk=PlyMatK{int}(k,:);
     vS_y{k}(rSk)=vS_sol{k};     % extension to (y,x_N\S).
     rcpq{k}=abs(vS_y{k}-x)<tol;
   elseif strcmp(str,'HMS_PN')
     vk=cell2mat(vS{int}(k));
     if length(vk)==1
        vS_sol{k}=PreKernel(vk); % solution y restricted to S.
        rSk=PlyMatK{int}(k,:);
        vS_y{k}(rSk)=vS_sol{k}; % extension to (y,x_N\S).
        rcpq{k}=abs(vS_y{k}-x)<tol;
      else
        try
          vS_sol{k}=Prenucl(vk); % solution y restricted to S.
        catch
          vS_sol{k}=PreNucl(vk); % use a third party solver instead!
        end
        rSk=PlyMatK{int}(k,:);
        vS_y{k}(rSk)=vS_sol{k}; % extension to (y,x_N\S).
        rcpq{k}=abs(vS_y{k}-x)<tol;
      end
   else % default
     vk=cell2mat(vS{int}(k));
     vS_sol{k}=PreKernel(vk); % solution y restricted to S.
     rSk=PlyMatK{int}(k,:);
     vS_y{k}(rSk)=vS_sol{k}; % extension to (y,x_N\S)
     rcpq{k}=PrekernelQ(clv,vS_y{k});
   end
  rcpQ{1,int}(k)=all(rcpq{1,k});
 end
KvS_sol{int,:}=vS_sol;
KvS_y{int,:}=vS_y;
krcpQ(int)=all(rcpQ{1,int});
end

kRCPQ=all(all(krcpQ));
%Formatting Output
if nargout>1
 kRCP=struct('kRCPQ',kRCPQ,'krcpQ',krcpQ);
 kRCPC={'vS',vS,'KvS_sol',KvS_sol,'KvS_y',KvS_y};
else
  kRCP=struct('kRCPQ',kRCPQ,'krcpQ',krcpQ);    
end

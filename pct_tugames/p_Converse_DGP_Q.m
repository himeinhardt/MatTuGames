function [CDGP CDGPC]=p_Converse_DGP_Q(v,x,str,tol)
% P_CONVERSE_DGP_Q checks whether an imputation x satisfies the converse derived game property 
% (CDGP) using Matlab's PCT, that is, a modified converse reduced game property (CRGP).
%
% Source: Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage [CDGP CDGPC]=p_Converse_DGP(v,x,str,tol)
% Define variables:
%  output: Fields
%  Q        -- Returns 1 (true) whenever the CDGP is satisfied, 
%              otherwise 0 (false).
%  cdgpQ    -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the derived game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell derived games on S 
%              with two players at x.
%  sV_x     -- Returns a vector of extended solutions x_s to x_N for 
%              all reduced games vS.              
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
  x=ModPreKernel(v);
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
PlyMat=false(N,n);
parfor k=1:n, PlyMat(:,k) = bitget(S,k)==1;end

% Checking now whether the imputation solves the derived game for 
% every two player coalitions.
sumPM=PlyMat*ones(n,1);
slcl2=sumPM==2;
PlyMat2=PlyMat(slcl2,:);
siPM2=size(PlyMat2,1);
Jmat=zeros(siPM2,2);

J=1:n;
for k=1:siPM2 
   Jmat(k,:)=J(PlyMat2(k,:));
end
Jmat=Jmat-1;
pw=2.^Jmat;
cl2=(pw*ones(2,1))';


vS=cell(1,siPM2);
stdsol=cell(1,siPM2);
crgpq=cell(1,siPM2);
crgpQ=false(1,siPM2);
sV_x=cell(1,siPM2);
rS=cell(1,siPM2);


parfor k=1:siPM2
sV_x{k}=x;
 if strcmp(str,'SHAP')
   vS{k}=HMS_Derived_game(v,x,cl2(k)); %Hart-MasColell derived game.
   stdsol{k}=StandardSolution(vS{k}); % solution x restricted to S.
   rS{k}=PlyMat2(k,:);
   sV_x{k}(rS{k})=stdsol{k}; % extension to (x,x_N\S).
   crgpq{k}=abs(sV_x{k}-x)<tol;
   crgpQ(k)=all(crgpq{k});
 else
   vS{k}=DerivedGame(v,x,cl2(k)); % Davis-Maschler derived game.
   stdsol{k}=StandardSolution(vS{k}); % solution x restricted to S.
   rS{k}=PlyMat2(k,:);
   sV_x{k}(rS{k})=stdsol{k}; % extension to (x,x_N\S).
   if strcmp(str,'MPRK')
     crgpQ(k)=ModPrekernelQ(v,sV_x{k});
   elseif strcmp(str,'PMPRK')
     crgpQ(k)=PModPrekernelQ(v,sV_x{k});
   elseif strcmp(str,'MODIC')
     crgpQ(k)=modiclusQ(v,sV_x{k});
   elseif strcmp(str,'PRK')
     crgpQ(k)=PrekernelQ(v,sV_x{k});
   elseif strcmp(str,'PRN')
    crgpq{k}=abs(sV_x{k}-x)<tol;
    crgpQ(k)=all(crgpq{k});
   else
    crgpQ(k)=ModPrekernelQ(v,sV_x{k});
   end
 end
end
CrgpQ=all(crgpQ);
%Formatting Output 
if nargout>1
 CDGP=struct('Q',CrgpQ,'cdgpQ',crgpQ);
 CDGPC={'vS',vS,'sV_x',sV_x};
else
  CDGP=struct('Q',CrgpQ,'cdgpQ',crgpQ);
end

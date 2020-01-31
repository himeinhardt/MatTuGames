function [COSED COSEDC]=SedcoconsQ(v,x,str,tol)
% SEDCOCONSQ checks whether an imputation x satisfies small excess difference converse 
% consistency (SEDCOCONS), that is, the converse SEDCONS (converse modified consistency).
%
% Source:  H. I. Meinhardt. Reconsidering Related Solutions of the Modiclus. Technical report, Karlsruhe Institute of Technology (KIT),
%                Karlsruhe, Germany, 2018c. URL http://dx.doi.org/10.13140/RG.2.2.27739.82729.
%          H. I. Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage [COSED COSEDC]=SedcoconsQ(v,x,str,tol)
%
% Define variables:
%
%  output: Fields
%  CrgpQ    -- Returns 1 (true) whenever SEDCOCONS is satisfied, 
%              otherwise 0 (false).
%  crgpQ    -- Gives a precise list of reduced games for which the 
%              restriction of x on S is a solution of the reduced game vS. 
%              It returns a list of zeros and ones.
%  vS       -- All Davis-Maschler or Hart-MasColell reduced games on S 
%              with two players at x.
%  sV_x     -- Returns a vector of extended solutions x_s to x_N for 
%              all reduced games vS.              
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'MAPRK' that is, the Davis-Maschler anti-reduced game 
%               in accordance with the modified anti-pre-kernel solution.
%              'PMAPRK' that is, the Davis-Maschler anti-reduced game 
%               in accordance with the proper modified anti-pre-kernel solution.
%              'SHAP' that is, Hart-MasColell anti-reduced game 
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler anti-reduced game
%               equivalence in accordance with the modiclus.
%              'APRN' that is, the Davis-Maschler anti-reduced game 
%               in accordance with the anti-pre-nucleolus.
%              'APRK' that is, the Davis-Maschler anti-reduced game 
%               in accordance with anti-pre-kernel solution.
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
%   03/10/2018        1.0             hme
%                


if nargin<2
  x=Anti_ModPreKernel(v);
  n=length(x);
  tol=10^6*eps;
  str='MAPRK';
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='MAPRK';
elseif nargin<4
  n=length(x);
  tol=10^6*eps;
else
  n=length(x);
end

N=length(v);
S=1:N;
PlyMat=false(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k)==1;end

% Checking now whether the imputation solves the reduced game for 
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

v_x=ECFloorGame(v,x);

vS=cell(1,siPM2);
stdsol=cell(1,siPM2);
crgpq=cell(1,siPM2);
crgpQ=false(1,siPM2);
sV_x=cell(1,siPM2);
rS=cell(1,siPM2);


for k=1:siPM2
sV_x{1,k}=x;
 if strcmp(str,'SHAP')
   vS{1,k}=HMS_RedGame(v_x,x,cl2(k)); %Hart-MasColell reduced game.
   stdsol{1,k}=StandardSolution(vS{1,k}); % solution x restricted to S.
   rS{k}=PlyMat2(k,:);
   sV_x{1,k}(rS{k})=stdsol{1,k}; % extension to (x,x_N\S).
   crgpq{k}=abs(sV_x{1,k}-x)<tol;
   crgpQ(k)=all(crgpq{k});
 else
   vS{1,k}=RedGame(v_x,x,cl2(k)); % Davis-Maschler reduced game.
   stdsol{1,k}=StandardSolution(vS{1,k}); % solution x restricted to S.
   rS{k}=PlyMat2(k,:);
   sV_x{1,k}(rS{k})=stdsol{1,k}; % extension to (x,x_N\S).
   if strcmp(str,'MAPRK')
     crgpQ(k)=Anti_ModPrekernelQ(v,sV_x{1,k});
   elseif strcmp(str,'PMAPRK')
     crgpQ(k)=Anti_PModPrekernelQ(v,sV_x{1,k});
   elseif strcmp(str,'MODIC')
     crgpQ(k)=modiclusQ(v,sV_x{1,k});
   elseif strcmp(str,'APRK')
     crgpQ(k)=Anti_PrekernelQ(v,sV_x{1,k});
   elseif strcmp(str,'APRN')
    crgpq{k}=abs(sV_x{1,k}-x)<tol;
    crgpQ(k)=all(crgpq{k});
   else
    crgpQ(k)=Anti_ModPrekernelQ(v,sV_x{1,k});
   end
 end
end
CrgpQ=all(crgpQ);
%Formatting Output 
if nargout>1
 COSED=struct('SEDCOCQ',CrgpQ,'sedcocQ',crgpQ);
 COSEDC={'vS',vS,'sV_x',sV_x};
else
  COSED=struct('SEDCOCQ',CrgpQ,'sedcocQ',crgpQ);
end

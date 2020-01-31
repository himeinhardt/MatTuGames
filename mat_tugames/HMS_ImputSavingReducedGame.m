function v_t=HMS_ImputSavingReducedGame(v,x,str)
% HMS_IMPUTSAVINGREDUCEDGAME computes from (v,x) all Hart/Mas-Colell imputation saving reduced 
% games on S at x of game v.
%
% Usage: v_t=HMS_ImputSavingReducedGame(v,x,str)
% Define variables:
%  output:
%  v_t{1,:} -- All Hart-MasColell imputation saving reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding Shapley values of all imputation saving reduced games.
%  v_t{3,:} -- The corresponding sub-coalitions which define a imputation saving reduced game.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'NUC' that is, the Hart-MasColell imputation saving reduced game 
%               in accordance with the nucleolus to compute 
%               a full set of reduced games.
%              'KRN' that is, the Hart-MasColell imputation saving reduced game 
%               in accordance with kernel solution to compute
%               a full set of reduced games.
%              'SHAP' that is, Hart-MasColell imputation saving reduced game 
%               in accordance with the Shapley Value, 
%               which is, the original definition. 
%              Default is 'SHAP'.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/13/2015        0.7             hme
%   02/23/2017        0.9             hme
%                

if nargin<2
  x=ShapleyValue(v);
  n=length(x);
  str='SHAP';
elseif nargin<3
  n=length(x);
  str='SHAP';
elseif nargin<4
  n=length(x);
else
  n=length(x);
end


N=length(v);  
v_t=cell(3,N-1);

for S=1:N-1
  [v_t{1,S} v_t{2,S} v_t{3,S}]=hms_impsavred_game(v,x,S,n,str);
end

%---------------------------------
function [vt subg_sh T]=hms_impsavred_game(v,x,S,n,str)


J=1:n;
plS=bitget(S,J)==1;
lmcS=bitget(S,J)==0;
plcS=J(lmcS);
cSpot=2.^(plcS-1);
cS=cSpot*ones(length(plcS),1);
T=SubSets(S,n);
lgt=length(T);
 

vt=zeros(1,lgt);
TorcS=bitor(T,cS);

subT=cell(1,lgt);
subg=cell(1,lgt);
subg_sh=cell(1,lgt);
lg=cell(1,lgt);
plT=cell(1,lgt);
Tz=cell(1,lgt);
sum_py=cell(1,lgt);

for k=1:lgt
 subT{k}=SubSets(TorcS(k),n);
 subg{k}=v(subT{k});
 if strcmp(str,'SHAP')
   subg_sh{k}=ShapleyValue(subg{k});
 elseif strcmp(str,'NUC')
   if length(subg{k})==1
      subg_sh{k}=subg{k};
   else
      try
        subg_sh{k}=cplex_nucl(subg{k});
      catch
        subg_sh{k}=nucl(subg{k}); % use a thrid party solver instead! 
      end
   end
 elseif strcmp(str,'KRN')
   subg_sh{k}=Kernel(subg{k});
 else
   subg_sh{k}=ShapleyValue(subg{k});
 end
 Qk=TorcS(k);
 it=0:-1:1-n;
 lg{k}=rem(floor(Qk(:)*pow2(it)),2)==1;
 plT{k}=J(lg{k});
 Tz{k}=ismember(plT{k},plcS);
 sum_py{k}=Tz{k}*subg_sh{k}';
 vt(k)=v(TorcS(k))-sum_py{k};
end

%% Constructing HMS Imputation Saving Reduced Game
sli=J(plS);
iS=2.^(sli-1);
for k=1:lgt
    cTi=T(k)==iS;
    if any(cTi)
       pl=sli(cTi);
       vxt(k)=min(x(pl),vt(k)); 
    else    
       vxt(k)=vt(k);
    end   
end

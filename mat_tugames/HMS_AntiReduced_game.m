function v_t=HMS_AntiReduced_game(v,x,str)
% HMS_ATNIREDUCED_GAME computes from (v,x) all Hart/Mas-Colell anti-reduced 
% games on S at x of game v.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: v_t=HMS_AntiReduced_game(v,x,str)
%
% Define variables:
%  output:
%  v_t{1,:} -- All Hart-MasColell anti-reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding Shapley values of all anti-reduced games.
%  v_t{3,:} -- The corresponding sub-coalitions which define an anti-reduced game.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'APRN' that is, the Hart-MasColell anti-reduced game 
%               in accordance with the pre-nucleolus.
%              'APRK' that is, the Hart-MasColell anti-reduced game 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell anti-reduced game 
%               in accordance with the Shapley Value, 
%               which is, the original definition. 
%              'MODIC' that is, the Hart-MasColell anti-reduced game 
%               in accordance with the modiclus.
%              Default is 'SHAP'.

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
  [v_t{1,S} v_t{2,S} v_t{3,S}]=hms_AntiRed_game(v,S,n,str);
end

%---------------------------------
function [vt subg_sh T]=hms_AntiRed_game(v,S,n,str)


J=1:n;
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
 elseif strcmp(str,'APRN')
   if length(subg{k})==1
      subg_sh{k}=subg{k};
   else
      try
        subg_sh{k}=cplex_AntiPreNucl_llp(subg{k});
      catch
        subg_sh{k}=Anti_PreNucl_llp(subg{k}); % use a thrid party solver instead! 
      end
   end
elseif strcmp(str,'MODIC')
   if length(subg{k})==1
      subg_sh{k}=subg{k};
   else
      try
        subg_sh{k}=cplex_modiclus(subg{k});
      catch
        subg_sh{k}=Modiclus(subg{k}); % use a thrid party solver instead! 
      end
   end
 elseif strcmp(str,'APRK')
   subg_sh{k}=Anti_PreKernel(subg{k});
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

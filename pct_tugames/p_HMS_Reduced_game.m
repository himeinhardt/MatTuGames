function v_t=p_HMS_Reduced_game(v,x,str)
% P_HMS_REDUCED_GAME computes from (v,x) all Hart/Mas-Colell reduced 
% games on S at x of game v. Using Matlab's PCT.
%
% Usage: v_t=p_HMS_Reduced_game(v,x,str)
% Define variables:
%  output:
%  v_t{:,1} -- All Hart-MasColell reduced games w.r.t. x.
%  v_t{:,2} -- The corresponding Shapley values of all reduced games.
%  v_t{:,3} -- The corresponding sub-coalitions which define a reduced game.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the Hart-MasColell reduced game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Hart-MasColell reduced game 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game 
%               in accordance with the Shapley Value, 
%               which is, the original definition. 
%              'MODIC' that is, the Hart-MasColell reduced game 
%               in accordance with the modiclus.
%              Default is 'SHAP'.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/21/2011        0.1 alpha        hme
%   06/17/2012        0.2 beta         hme
%   10/27/2012        0.3              hme
%   05/16/2014        0.5              hme
%   02/10/2018        0.9              hme
%   04/01/2020        1.9              hme
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
v1_t=cell(1,N-1);
v2_t=cell(1,N-1);
v3_t=cell(1,N-1);

parfor S=1:N-1
  [v1_t{1,S} v2_t{1,S} v3_t{1,S}]=hms_red_game(v,S,n,str);
end


v_t={v1_t v2_t v3_t};
%v_t=[v1_t v2_t v3_t];

%---------------------------------
function [vt subg_sh T]=hms_red_game(v,S,n,str)

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
 elseif strcmp(str,'PRN')
   if length(subg{k})==1
      subg_sh{k}=subg{k};
   else
      try
        subg_sh{k}=cplex_prenucl_mod4(subg{k});
      catch
        subg_sh{k}=PreNucl(subg{k});
      end
   end
 elseif strcmp(str,'PRK')
   subg_sh{k}=PreKernel(subg{k});
 elseif strcmp(str,'MODIC')
   if length(subg{k})==1
      subg_sh{k}=subg{k};
   else
      try
        subg_sh{k}=cplex_modiclus(subg{k});
      catch
        subg_sh{k}=Modiclus(subg{k});
      end
   end
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

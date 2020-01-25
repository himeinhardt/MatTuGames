function v_xt=ImputSavingReducedGame(v,x)
% ImputSavingReducedGame computes from (v,x) all imputation saving reduced games on S at x of
% game v.
%
% Usage: v_t=ImputSavingReducedGame(v,x) 
% Define variables:
%  output:
%  v_xt{1,:} -- All imputation saving reduced games w.r.t. x.
%  v_xt{2,:} -- All Davis-Maschler reduced games w.r.t. x.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/12/2015        0.7             hme
%                


N=length(v);  
n=length(x);
v_xt=cell(2,N-1);

for S=1:N-1
  [v_xt{1,S} v_xt{2,S}]=impsav_game(v,x,S,n);
end

%---------------------------------
function [vxt vt]=impsav_game(v,x,S,n)

J=1:n;
plS=bitget(S,J)==1;
lmcS=plS==0;
plcS=J(lmcS);
cSpot=2.^(plcS-1);
cS=cSpot*ones(length(plcS),1);
T=SubSets(S,n);
lgt=length(T);
vt=zeros(1,lgt);
vxt=zeros(1,lgt);
y=x(lmcS);
Q=SubSets(cS,n);
PayQ=additive_game(y);
TorQ=cell(lgt-1);
%% Constructing DM Reduced Game
for k=1:lgt-1;
   TorQ{k}=bitor(T(k),Q);
   vt(k)=max(v(TorQ{k})-PayQ);
   vt(k)=max(v(T(k)),vt(k)); % Considering the empty set of Q.
end
vt(lgt)=x*plS';

%% Constructing Imputation Saving Reduced Game
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

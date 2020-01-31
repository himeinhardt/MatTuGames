function v_t=p_DM_AntiReduced_game(v,x)
% P_DM_ANTIREDUCED_GAME computes from (v,x) all anti-reduced games on S at x of
% game v using Matlab's PCT.
%
% Usage: v_t=p_DM_AntiReduced_game(v,x) 
% Define variables:
%  output:
%  v_t{1,:} -- All Davis-Maschler anti-reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding sub-coalitions which define an anti-reduced game.
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
%-----------------------------------------------
%   03/01/2018        1.0             hme
%                


N=length(v);  
n=length(x);
v1_t=cell(1,N-1);
v2_t=cell(1,N-1);

spmd
 codistributed(v1_t);
 codistributed(v2_t);
end

parfor S=1:N-1
  [v1_t{1,S} v2_t{1,S}]=AntiRed_game(v,x,S,n);
end
v_t={v1_t, v2_t};

%---------------------------------
function [vt T]=AntiRed_game(v,x,S,n)

J=1:n;
plS=bitget(S,J);
lmcS=plS==0;
plcS=J(lmcS);
cSpot=2.^(plcS-1);
cS=cSpot*ones(length(plcS),1);
T=SubSets(S,n);
lgt=length(T);
vt=zeros(1,lgt);
Q=SubSets(cS,n);
it=0:-1:1-n;
plmQ=rem(floor(Q(:)*pow2(it)),2);
PayQ=plmQ*x';

TorQ=cell(lgt-1);

for k=1:lgt-1;
   TorQ{k}=bitor(T(k),Q);
   vt(k)=min(v(TorQ{k})-PayQ');
   vt(k)=min(v(T(k)),vt(k)); % Considering the empty set of Q.
end
vt(end)=x*plS';

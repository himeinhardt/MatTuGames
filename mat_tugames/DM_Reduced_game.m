function v_t=DM_Reduced_game(v,x)
% DM_REDUCED_GAME computes from (v,x) all reduced games on S at x of
% game v.
%
% Usage: v_t=DM_Reduced_game(v,x) 
% Define variables:
%  output:
%  v_t{1,:} -- All Davis-Maschler reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding sub-coalitions which define a reduced game.
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
%   08/19/2010        0.1 beta        hme
%   06/27/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%                


N=length(v);  
n=length(x);
v_t=cell(2,N-1);

for S=1:N-1
  [v_t{1,S} v_t{2,S}]=red_game(v,x,S,n);
end

%---------------------------------
function [vt T]=red_game(v,x,S,n)

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
   vt(k)=max(v(TorQ{k})-PayQ');
   vt(k)=max(v(T(k)),vt(k)); % Considering the empty set of Q.
end
vt(end)=x*plS';

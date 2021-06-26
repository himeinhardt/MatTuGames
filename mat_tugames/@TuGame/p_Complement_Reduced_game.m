function v_t=p_Complement_Reduced_game(clv,x)
% P_COMPLEMENT_REDUCED_GAME computes from (v,x) all complement reduced games on S at x of
% game v using Matlab's PCT.
%
% Source: Moulin H (1985) The separability axiom and equal sharing methods. J of Econ Theory 36:120-148
%
%
% Usage: v_t=clv.p_Complement_Reduced_game(x) 
% Define variables:
%  output:
%  v_t{1,:} -- All complement reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding sub-coalitions which define a complement reduced game.
%
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/13/2020        1.9             hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
v1_t=cell(1,N-1);
v2_t=cell(1,N-1);

spmd
 codistributed(v1_t);
 codistributed(v2_t);
end

parfor S=1:N-1
  [v1_t{1,S} v2_t{1,S}]=cmp_red_game(v,x,S,N,n);
end

v_t={v1_t, v2_t};

%---------------------------------
function [vt T]=cmp_red_game(v,x,S,N,n)

CMP=N-S;	
J=1:n;
plC=bitget(CMP,J);
T=SubSets(S,n);
lgt=length(T);
vt=zeros(1,lgt);
PayC=plC*x';
TorQ=cell(lgt);

for k=1:lgt;
   TorC=bitor(T(k),CMP);
   vt(k)=v(TorC)-PayC';
end

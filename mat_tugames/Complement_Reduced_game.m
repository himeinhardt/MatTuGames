function v_t=Complement_Reduced_game(v,x)
% COMPLEMENT_REDUCED_GAME computes from (v,x) all complement reduced games on S at x of
% game v.
%
% Source: Moulin H (1985) The separability axiom and equal sharing methods. J of Econ Theory 36:120-148
%
%
% Usage: v_t=Complement_Reduced_game(v,x) 
% Define variables:
%  output:
%  v_t{1,:} -- All complement reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding sub-coalitions which define a complement reduced game.
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
%   06/13/2020        1.9             hme
%                


N=length(v);  
n=length(x);
v_t=cell(2,N-1);

for S=1:N-1
  [v_t{1,S} v_t{2,S}]=cmp_red_game(v,x,N,S,n);
end

%---------------------------------
function [vt T]=cmp_red_game(v,x,N,S,n)

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

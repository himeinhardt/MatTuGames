function [vS,T]=CmpRedGame(clv,x,S)
% CMPREDGAME computes from (v,x,S) a complement reduced game vS on S at x for
% game v.
%
% Source: Moulin H (1985) The separability axiom and equal sharing methods. J of Econ Theory 36:120-148
% 
%
% Usage: [vS T]=CmpRedGame(clv,x,S)
% Define variables:
%  output:
%  vS      -- The complement reduced game vS w.r.t. x.
%  T       -- The corresponding sub-coalitions of S which define 
%             the complement reduced game vS.
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n).
%  S        -- A coalition/set identified by its unique integer representation.
%


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
n=clv.tuplayers;
N=clv.tusize;
CMP=N-S;
J=1:n;
plC=bitget(CMP,J);
T=SubSets(S,n);
lgt=length(T);
vS=zeros(1,lgt);
PayC=plC*x';
TorQ=cell(lgt);

for k=1:lgt;
   TorC=bitor(T(k),CMP);
   vS(k)=v(TorC)-PayC';
end



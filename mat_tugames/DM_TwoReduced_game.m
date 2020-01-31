function [v_t,sS,PlyMat2]=DM_TwoReduced_game(v,x)
% DM_TWOREDUCED_GAME computes from (v,x) all single and two-person reduced games on S at x of
% game v.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: v_t=DM_TwoReduced_game(v,x) 
% Define variables:
%  output:
%  v_t{1,:} -- All Davis-Maschler single and two-person reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding sub-coalitions which define a reduced game.
%  sS       -- The set of singleton and two-person coalitions.
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
%   08/13/2018        1.0             hme
%                


N=length(v);  
n=length(x);
S=1:N;
PlyMat=false(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k)==1;end

sumPM=PlyMat*ones(n,1);
slcl2=sumPM<=2;
sS=S(slcl2);
PlyMat2=PlyMat(slcl2,:);
lS2=length(sS);

v_t=cell(2,lS2);

for k=1:lS2
  [v_t{1,k} v_t{2,k}]=red_game(v,x,sS(k),n);
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

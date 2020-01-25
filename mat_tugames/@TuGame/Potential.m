function [pot,sh]=Potential(clv)
% POTENTIAL computes the potential and the Shapley value of a TU-game v.
%
% Usage: [pot,sh]=Potential(clv)
%
% Define variables:
%  output:
%  pot      -- The potential of the game v.
%  sh       -- The Shapley value of a TU-game v.
%
%  input:
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/07/2014        0.5             hme
%

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
Si=clv.tuSi;

if N==1
  sh=v;return;
 else
end

% Converting logical to double to speed up computation
% for weighted majority games.
if isnumeric(v)==0
   v=double(v);
end

pot=zeros(1,N);
k=1:n;
sC=2.^(k-1);
pot(sC)=v(sC);
S=1:N;
tv=true(1,N);
tv(sC)=false;
cpl=S(tv);
trv=v(cpl);
lc=length(cpl);
for i=1:lc
  pot=potential2(trv(i),pot,cpl(i),k);
end

sh=pot(N)-pot(Si);


%----------------------------
function p=potential2(vS,p,S,J)

clS=bitget(S,J)==1;
plS=J(clS);
Sni=S-2.^(plS-1);
lgS=length(Sni);
p(S)=(vS+p(Sni)*ones(lgS,1))/lgS;

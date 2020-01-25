function [dsh,pot]=DiscShapleyValue(v,alp)
% DISCSHAPLEYVALUE computes the discounted Shapley value and potential of a TU-game v.
%
% Usage: [dsh,pot]=DiscShapleyValue(v,alp)
% Define variables:
%  output:
%  pot       -- The discounted potential of the game v.
%  dsh       -- The discounted Shapley value of a TU-game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%  alp      -- A discount value, a number between (0,1)
%              if alp=0, we get
%                 the equal division.
%              if alp=1, then we get
%                 the Shapley value.
%              Default value 1.   
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/31/2014        0.5             hme
%

if nargin <2
   alp=1;
end
N=length(v);
[~, n]=log2(N);
if N==1
  dsh=v;return;
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
  pot=potential2(trv(i),pot,cpl(i),k,alp);
end

Nk=bitset(N,k,0);
dsh=pot(N)-alp*pot(Nk);


%----------------------------
function p=potential2(vS,p,S,J,alp)

clS=bitget(S,J)==1;
plS=J(clS);
Sni=S-2.^(plS-1);
lgS=length(Sni);
p(S)=(vS+alp*p(Sni)*ones(lgS,1))/lgS;

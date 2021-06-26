function epsd=EPSDValue(v)
% EPSDValue computes the egalitarian proportional surplus division value 
% of a individually positive TU-game.
%
% Usage: epsd=EPSDValue(v)
%
% Define variables:
%  output:
%  pesd      -- Egalitarian proportional surplus division 
%               value of a individually positive TU-game. 
%
%  input:
%  v         -- A TU-Game of length 2^n-1.
%
%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/27/2020        1.9             hme
%
N=length(v);
[~, n]=log2(N);
if N==1
  esd=v;return;
 else
end
k=1:n;
sC=2.^(k-1);
vi=v(sC);
Svi=sum(vi);
ipQ=all(vi>0);
if ipQ==1
   epsd=(1/n)*Svi + (vi/Svi)*(v(N)-Svi);
else
   epsd=-inf(1,n); 
   warning("EPSDValue: Game is not individually positive!!");  
end    
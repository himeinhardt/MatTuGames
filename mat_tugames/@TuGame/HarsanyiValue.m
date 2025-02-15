function hv=HarsanyiValue(clv)
% HARSANYIVALUE computes a Harsanyi-value of a TU-game v.
%
% Usage: hv=clv.HarsanyiValue()
%
% Define variables:
%  output:
%  hv       -- A Harsanyi-value of a TU-game v.
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
%   01/28/2021        1.9.1           hme
%
N=clv.tusize;
n=clv.tuplayers;    
hd_v=clv.harsanyi_dividends();
S=1:N;
PlyMat=zeros(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k);end    
pmsz=sum(PlyMat,2);
probm=(PlyMat./pmsz)';
hv=(probm*hd_v')';

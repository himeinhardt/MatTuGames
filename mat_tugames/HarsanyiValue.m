function hv=HarsanyiValue(v)
% HARSANYIVALUE computes a Harsanyi-value of a TU-game v.
%
% Usage: hv=HarsanyiValue(v)
%
% Define variables:
%  output:
%  hv       -- A Harsanyi-value of a TU-game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/12/2021        1.9.1           hme
%
    
N=length(v);
[~, n]=log2(N);    
hd_v=harsanyi_dividends(v);
S=1:N;
PlyMat=zeros(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k);end    
pmsz=sum(PlyMat,2);
probm=(PlyMat./pmsz)';
hv=(probm*hd_v')';
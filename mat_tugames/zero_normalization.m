function zv=zero_normalization(v)
% ZERO_NORMALIZATION computes from the game v  
% the corresponding zero-normalized game zv.
%
% Usage: zv=zero_normalization(v)
% Define variables:
%  output:
%  zv       -- The zero-normalized of the TU-game v.
%  input:
%  v        -- A TU-game v of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/13/2010        0.1 beta        hme
%   10/27/2012        0.3             hme
%                

N=length(v);
[~, n]=log2(N);

k=1:n;
sC=bitset(0,k);
svi=v(sC);

sumvi=zeros(1,N);
sumvi=svi(1); for ii=2:n, sumvi=[sumvi svi(ii) sumvi+svi(ii)]; end

zv=(v-sumvi);

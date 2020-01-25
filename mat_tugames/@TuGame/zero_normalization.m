function zv=zero_normalization(clv)
% ZERO_NORMALIZATION computes from the game v  
% the corresponding zero-normalized game zv.
%
% Usage: zv=zero_normalization(clv)
%
% Define variables:
%  output:
%  zv       -- The zero-normalized of the TU-game v.
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
%   10/28/2012        0.3             hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
k=1:n;
sC=bitset(0,k);
svi=v(sC);

sumvi=zeros(1,N);
sumvi=svi(1); for ii=2:n, sumvi=[sumvi svi(ii) sumvi+svi(ii)]; end

zv=(v-sumvi);

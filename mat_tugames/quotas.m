function [qt, vni]=quotas(v)
%QUOTAS computes the quotas of a tu-game.
%
% Usage: qt=quotas(v)
%
% Define variables:
%  output:
%  qt       -- the quotas of a tu-game v.
%  vni      -- the coalitional values of n-1.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/10/2013        0.4             hme
%                



N=length(v);
[~, n]=log2(N);
k=1:n;
Si=bitset(N,k,0);
vni=v(Si);
it=0:-1:1-n;
mat=rem(floor(Si(:)*pow2(it)),2);
qt=(mat\vni')';

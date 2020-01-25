function [qt, vni]=quotas(clv)
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
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/10/2013        0.4             hme
%                


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
Si=clv.tuSi;

vni=v(Si);
it=0:-1:1-n;
mat=rem(floor(Si(:)*pow2(it)),2);
qt=(inv(mat)*vni')';

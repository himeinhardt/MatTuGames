function zov=ZeroOne_Normalization(v)
% ZEROONE_NORMALIZATION computes from the game v  
% the corresponding (0,1)-normalized game zov.
%
% Usage: zov=ZeroOne_Normalization(v)
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
%   10/28/2012        0.3             hme
%                

zv=zero_normalization(v);
if zv(end)~=0
   zov=zv/zv(end);
 else
   warning('ZN:zr','No (0,1)-normalized game computed, otherwise division by zero.');
end

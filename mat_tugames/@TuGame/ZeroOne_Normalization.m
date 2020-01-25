function zov=ZeroOne_Normalization(clv)
% ZEROONE_NORMALIZATION computes from the game v  
% the corresponding (0,1)-normalized game zov.
%
% Usage: zov=ZeroOne_Normalization(clv)
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

zv=zero_normalization(clv);
if zv(end)~=0
   zov=zv/zv(end);
 else
   warning('ZN:zr','No (0,1)-normalized game computed, otherwise division by zero.');
end

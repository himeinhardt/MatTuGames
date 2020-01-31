function ov=One_Normalization(v)
% ONE_NORMALIZATION computes from the game v  
% the corresponding one-normalized game ov.
%
% Usage: zov=ZeroOne_Normalization(v)
% Define variables:
%  output:
%  ov       -- The zero-normalized of the TU-game v.
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
%   02/28/2016        0.8             hme
%                

if v(end)~=0
   ov=v/v(end);
 else
   warning('ZN:zr','No 1-normalized game computed, otherwise division by zero.');
end

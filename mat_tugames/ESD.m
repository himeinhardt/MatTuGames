function esd=ESD(v)
% ESD computes the equal surplus division of a TU-game.
%
% Usage: esd=ESD(v)
%
% Define variables:
%  output:
%  esd      -- Equal surplus division of a TU-game. 
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/13/2014        0.5             hme
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
esd=vi+(v(N)-sum(vi))/n;

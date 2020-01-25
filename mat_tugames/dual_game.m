function dv=dual_game(v)
% DUAL_GAME computes from the game v its dual game dv.
%
% Usage: dv=dual_game(v)
% Define variables:
%  output:
%  dv       -- The dual game of v.
%  input:
%  v        -- A TU-game.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/06/2010        0.1 beta        hme
%                
N=length(v);
S=1:N-1;
CN=N-S;
cv=v(CN);
cv(N)=0;
dv=(v(N)-cv);

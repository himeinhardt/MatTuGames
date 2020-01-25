function dv=dual_game(clv)
% DUAL_GAME computes from the game v its dual game dv.
%
% Usage: dv=dual_game(clv)
%
% Define variables:
%  output:
%  dv       -- The dual game of v.
%
%  input:
%  clv        -- TuGame class object.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%   
v=clv.tuvalues;
N=clv.tusize;             
S=1:N-1;
CN=N-S;
cv=v(CN);
cv(N)=0;
dv=(v(N)-cv);

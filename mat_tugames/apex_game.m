function sv=apex_game(ii,n)
% APEX_GAME determines an apex game.
%
% Usage: sv=apex_game(ii,n)
%
% Define variables:
%  output:
%  sv        -- The apex game (simple game) of lenght 2^n-1.
%
%  input:
%  ii        -- The apex or main player, must be an integer 
%               between 1 and n.    
%  n         -- The number of player involved in the apex game.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/16/2013        0.4             hme
% 
N=2^n-1;
S=1:N;
Si=bitget(S,ii)==1;
lSi=2^(ii-1);
Si(lSi)=0;
sv =Si;
Ni=bitset(N,ii,0);
sv(Ni)=1;
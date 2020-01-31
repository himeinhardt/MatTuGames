function sv=jury_game(th,n,k)
% JURY_GAME computes from a quota and the number of jurors involved  
% the corresponding simple game.
%
% Usage: sv=simple_game(th,n)
% Example:
% Let th=10 and n=12;
%    sv=simple_game(10,12);
% computes the corresponding simple game.
%
%
% Define variables:
%  output:
%  sv       -- A simple game.
%  input:
%  th       -- Number of jurors required to convict a defendant.
%  n        -- The number of players involved (positive number).
%  k        -- Veto juror (judge/player)
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/01/2015        0.6             hme
%


N=2^n-1;
S=1:N;
it=0:-1:1-n;
sS=rem(floor(S(:)*pow2(it)),2);
szc=sS*ones(n,1);
if nargin == 3
  kS=bitget(S,k);
  sv=(szc>=th)';
  sv = sv & kS;    
else
  sv=(szc>=th)';
end

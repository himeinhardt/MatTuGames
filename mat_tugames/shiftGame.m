function vt=shiftGame(v,t);
% SHIFTGAME computes from the game v the t-shift game of v.
% 
% Usage: vt=shiftGame(v,t)
%

% Define variables:
%  output:
%  vt       -- Returns the characteristic function of the 
%              t-shfit game of v.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  t        -- epsilon value (cost) to form a proper subcoalition. 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/19/2013        0.4             hme
%                


vt=v+t;
vt(end)=v(end);



function vt=shiftGame(clv,t);
% SHIFTGAME computes from the game v the t-shift game of v.
% 
% Usage: vt=clv.shiftGame(t)
%

% Define variables:
%  output:
%  vt       -- Returns the characteristic function of the 
%              t-shfit game of v.
%
%  input:
%  clv      -- TuGame class object.
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

v=clv.tuvalues;
N=clv.tusize;

vt=v+t;
vt(N)=v(N);



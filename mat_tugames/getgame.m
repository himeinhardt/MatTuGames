function v=getgame(hd,n);
% GETGAME computes a game from the unanimity coordinates hd 
%
% Usage: v=getgame(hd,n)
% Define variables:
%  output:
%  v        -- A TU-bankruptcy game.
%  input:
%  hd        -- The harsanyi dividends.
%  n         -- The number of player involved (optional)
%
%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/14/2011        0.1 beta        hme
%   10/27/2012        0.3             hme
%
%

if nargin<1
    error('At least the game must be given!');
elseif nargin<2
    N=length(hd);
    [~, n]=log2(N);
    if (2^n-1)~=N
      error('Unanimity coordinates have not the correct size!');
    end
else
end

gb=game_basis(n);
gb=sparse(gb);
v=gb*hd';
v=v';

function v=min_game(lam)
% MIN_GAME(lam) generates a minimum game. Every market game can be
% represented as a minimum game.
%
% Source: E. Kalai and E. Zemel, Generalized Network problems yielding
% totally balanced games, Operations Research Vol. 30, pp. 998-1008 (1982)).
% 
% Usage: v=min_game(lam)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  lam      -- A matrix of measures of size (r x n). 
%
%
% Example:
% Define a matrix of measures given by
% >> lam=[3 0 0 0;0 3 0 0;0 0 2 1]
%
% lam =
%
%     3     0     0     0
%     0     3     0     0
%     0     0     2     1
%
% Then the minimum game is specified by 
% >> v=min_game(lam)
%
% v =
%
%     0     0     0     0     0     0     2     0     0     0     1     0     0     0     3
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/19/2017        0.9             hme
%


[r,n]=size(lam);
N=2^n-1;
S=1:N;
it=0:-1:1-n;
mat=(rem(floor(S(:)*pow2(it)),2)==1);
ms=mat*lam';
for k=1:N
  v(k)=min(ms(k,:));
end  
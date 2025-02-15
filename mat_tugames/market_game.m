function v=market_game(P,Q,n,slc)
% MARKET_GAME generates a short-sided market game from a trading situation.  
% 
% Usage: v=market_game(P,Q,n)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  P        -- An integer to specify the set of producers. 
%  Q        -- An integer to specify the set of buyers s.t. P and Q
%              partition N.
%  n        -- An integer to specify the number of traders involved in
%              the game.
%  slc      -- Holding of units of the numeraire good per buyer, default is 1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/19/2013        0.4             hme
%   06/05/2023        1.9.1           hme
%
msg=nargchk(3,4,nargin);
error(msg);

if nargin < 4
    slc=1;
end    
N=2^n-1;
S=1:N;
if isempty(bitand(P,Q))==1
   error('Both sets must be disjoint!')
else
   if bitor(P,Q)~=N
     error('Both sets must partition the grand coaliton!')
   end
end 
msp=false(N,n);
msq=false(N,n);
v=zeros(1,N);
sp=bitand(S,P);
sq=bitand(S,Q);
for ii=1:n
       msp(:,ii)=bitget(sp,ii);
       msq(:,ii)=bitget(sq,ii);
end
psz=msp*ones(n,1);
qsz=slc*msq*ones(n,1);
v=min(psz,qsz)';




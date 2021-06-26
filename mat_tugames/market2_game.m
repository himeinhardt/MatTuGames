function v=market2_game(P,Q,n,csz,vsz,dfv,vN)
% MARKET2_GAME(P,Q,n) generates a producer and buyer game.
% 
% Usage: v=market2_game(P,Q,n,scz,vsz,dfv,vN) 
%
% Example:
% Set P=[1 3 5 7];
% Set Q=[2 4 6 8]
% n=8;
% csz=2;
% vsz=4;
% dfv=0;
% vN=16;
%
% Execute 
%     v=market2_game(P,Q,n,scz,vsz,dfv,vN);
% to get an eight person game
% or simply
%     v=market2_game;
% to get the default game.
%
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  P        -- A set of integers to specify the set of producers. 
%  Q        -- A set of integers to specify the set of buyers s.t. P and Q
%              partition N.
%  n        -- An integer to specify the number of persons involved in
%              the game.
%  csz      -- The size of buyers coalitions that join a producer.
%  vsz      -- Value that should assigned to such a coalition from csz.
%  dfv      -- default value for the remaining coalitions
%  vN       -- value of the grand coalition; 
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/22/2014        0.6             hme
%

if nargin < 1
   P=1:4;
   Q=5:8;
   n=8;
   csz=3;
   vsz=0;
   dfv=-4/7;
   vN=4;
elseif nargin<2
   Q=5:8;
   n=8;
   csz=3;
   vsz=0;
   dfv=-4/7;
   vN=4;
elseif nargin<3
   n=8;
   csz=3;
   vsz=0;
   dfv=-4/7;
   vN=4;
elseif nargin<4
   csz=3;
   vsz=0;
   dfv=-4/7;
   vN=4;
elseif nargin<5
   vsz=0;
   dfv=-4/7;
   vN=4;
elseif nargin<6
   dfv=-4/7;
   vN=4;
elseif nargin<7
   vN=4;
end

k=1:n;
vi=2.^(k-1);
lp=length(P);
N=2^n-1;
v=ones(1,N)*dfv;
v(vi)=vsz;
sC=SubCoalitions(Q,n,csz);
for k=1:lp
    v(bitset(sC,k))=vsz;
end
v(N)=vN;

function v=cardinality_game(n,k)
% CARDINALITY_GAME(n,k) assigns a zero to a coalition of size<=k<n,
% otherwise its cardinality.
% Usage: v=cardinality_game(n,k)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  n        -- An integer to specify the number of persons involved in
%              the game.
%  k        -- An integer to specify the coalition which gets zero. 
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/20/2011        0.1 beta        hme
%   10/18/2012        0.3             hme
%                
N=2^n-1;
S=1:N;
csz=zeros(1,N);
v=zeros(1,N);
PlyMat=false(N,1);
for ii=1:n, PlyMat(:,ii) = bitget(S,ii)==1;end
csz=PlyMat*ones(n,1);
lg1=csz>k;
cs=S(lg1);
v(cs)=csz(lg1);

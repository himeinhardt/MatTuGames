function v=production_game(n,k)
% production_GAME(n,k) assigns a zero to a coalition of size<=k<n,
% otherwise its cardinality-1.
% Usage: v=production_game(n,k)
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
%   07/20/2013        0.4             hme
%                
N=2^n-1;
S=1:N;
csz=zeros(1,N);
v=zeros(1,N);
for ii=1:n, PlyMat(:,ii) = bitget(S,ii);end
csz=PlyMat*ones(n,1);
lg1=csz>=k;
cs=S(lg1);
v(cs)=csz(lg1)-1;




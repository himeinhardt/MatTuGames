function v=oddeven_game(n)
% ODDEVEN_GAME(n) assigns |S|-1 if S is odd and 
% |S|+1 if S is even. 
% otherwise its cardinality.
% Usage: v=oddeven_game(n)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  n        -- An integer to specify the number of persons involved in
%              the game.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/16/2017        0.9             hme
%                
N=2^n-1;
S=1:N;
csz=zeros(1,N);
v=zeros(1,N);
for ii=1:n, PlyMat(:,ii) = bitget(S,ii);end
csz=PlyMat*ones(n,1);
cs=csz';
c1=cs(mod(cs,2)==0);
c2=cs(mod(cs,2)==1);
ev=S(mod(cs,2)==0);
odv=S(mod(cs,2)==1);
v(ev)=c1+1;
v(odv)=c2-1;





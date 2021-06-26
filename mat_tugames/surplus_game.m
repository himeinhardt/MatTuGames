function v=surplus_game(c)
%SURPLUS_GAME computes from a cost game c the corresponding surplus game v.
%
% Usage: v=surplus_game(c)
% Define variables:
%  output:
%  v        -- The surplus game v of length 2^n-1.
%  input:
%  c        -- A cost Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/11/2021        1.9             hme
%                



N=length(c);
[~, n]=log2(N);
v=zeros(1,N);
for S=1:N;
    v=surplus(c,v,S,n);
end

%--------------------------------
function v=surplus(c,v,S,n)

it=0:-1:1-n;
vecS=rem(floor(S(:)*pow2(it)),2)==1;;

J=1:n;
sP=J(vecS);
sC=bitset(0,sP);

lgS=length(sP);

v(S)=c(S)-(c(sC)*ones(lgS,1));

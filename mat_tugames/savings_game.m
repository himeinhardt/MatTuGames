function v=savings_game(c)
%  SAVINGS_GAME computes from a cost game c the corresponding savings game v.
%
% Usage: v=savings_game(c)
% Define variables:
%  output:
%  v        -- The savings game v of length 2^n-1.
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
%   08/12/2010        0.1 beta        hme
%   06/24/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%                



N=length(c);
[~, n]=log2(N);
v=zeros(1,N);

for S=1:N;
v=savings(c,v,S,n);
end

%--------------------------------
function v=savings(c,v,S,n)

it=0:-1:1-n;
vecS=rem(floor(S(:)*pow2(it)),2)==1;;

J=1:n;
sP=J(vecS);
sC=bitset(0,sP);

lgS=length(sP);

v(S)=(c(sC)*ones(lgS,1))-c(S);

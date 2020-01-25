function [vk kcQ]=p_k_cover(clv,k)
% P_K_COVER computes from the game v and the integer k the corresponding 
% k-game.
%
% Usage: [vk kcQ]=p_k_cover(clv,k)
%
% Define variables:
%  output:
%  vk       -- The derived k-game from game v. 
%  kcQ      -- Returns 1 (true) whenever vk is the k-cover of game v.
%
%  input:
%  clv      -- TuGame class object.
%  k        -- An integer, not greater than n.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/30/2012        0.3              hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

vk=zeros(1,N);
gv=zeros(1,N);
kcq=zeros(1,N);
bv=zeros(1,n);
lv=zeros(1,n);

[gv, bv, lv]=p_Gap(v);


for S=1:N;
 vk=vk_game(v,vk,gv,bv,S,k,n,N);
end
kcq=vk>=v;
kcQ=all(kcq);


%----------------------------
function [vk]=vk_game(v,vk,gv,bv,S,k,n,N)

it=0:-1:1-n;
vecS=rem(floor(S(:)*pow2(it)),2)==1;

J=1:n;
sP=J(vecS);


if length(sP)<k
    vk(S)=v(S);
 else
    lgS=length(sP);
    vk(S)=bv(sP)*ones(lgS,1)-gv(N);
end

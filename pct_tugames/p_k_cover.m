function [vk kcQ]=p_k_cover(v,k)
% P_K_COVER computes from the game v and the integer k the corresponding 
% k-game.
%
% Usage: [vk kcQ]=p_k_cover(v,k)
% Define variables:
%  output:
%  vk       -- The derived k-game from game v. 
%  kcQ      -- Returns 1 (true) whenever vk is the k-cover of game v.
%  input:
%  v        -- A TU-game of length 2^n-1.
%  k        -- An integer, not greater than n.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/20/2011        0.1 alpha        hme
%   06/19/2012        0.2 beta         hme
%   10/27/2012        0.3              hme
%                


N=length(v);
[~, n]=log2(N);
vk=zeros(1,N);
gv=zeros(1,N);
kcq=zeros(1,N);
bv=zeros(1,n);
lv=zeros(1,n);

[gv bv lv]=p_Gap(v);


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

function [g bv lv]=p_Anti_Gap(v)
% P_NTI_GAP computes the anti-gap function from game v using Matlab's PCT.
%
% Usage: [g bv lv]=p_Anti_Gap(v)
% Define variables:
%  output:
%  g        -- The anti-gap function of game v. A vector of length 2^n-1.
%  bv       -- The anti-upper vector/payoff of game v.
%  lv       -- The anti-lower/concession vector of game v.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/07/2021        1.9              hme
%                



N=length(v);
[~, n]=log2(N);
% anti-upper vector
Si=false(n,1);
k=1:n;
Si=bitset(N,k,0);
bv=v(N)-v(Si);
%S=1:N;

% Computing the anti-gap function w.r.t. v.
Bm=bv(1); for ii=2:n, Bm=[Bm bv(ii) Bm+bv(ii)]; end
g=v-Bm;

lv=zeros(1,n);

% anti-concession vector
S=1:N;
parfor i=1:n
   a_i=bitget(S,i)==1;
   lv(i)=min(g(a_i));
end

lv=gather(lv);

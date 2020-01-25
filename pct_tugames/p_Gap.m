function [g bv lv]=p_Gap(v)
% P_GAP computes the gap function from game v using Matlab's PCT.
%
% Usage: [g bv lv]=p_Gap(v)
% Define variables:
%  output:
%  g        -- The gap function of game v. A vector of length 2^n-1.
%  bv       -- The upper vector/payoff of game v.
%  lv       -- The lower/concession vector of game v.
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
%   05/21/2011        0.1 alpha        hme
%   09/14/2012        0.2              hme
%   10/27/2012        0.3              hme
%                



N=length(v);
[~, n]=log2(N);
% upper vector
Si=false(n,1);
k=1:n;
Si=bitset(N,k,0);
bv=v(N)-v(Si);
S=1:N;

% Computing the gap function w.r.t. v.
Bm=bv(1); for ii=2:n, Bm=[Bm bv(ii) Bm+bv(ii)]; end
g=Bm-v;

lv=zeros(1,n);

% concession vector
parfor i=1:n
   a_i=bitget(S,i)==1;
   lv(i)=min(g(a_i));
end

lv=gather(lv);

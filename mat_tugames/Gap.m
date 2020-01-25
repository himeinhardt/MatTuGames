function [g bv lv]=Gap(v)
% GAP computes the gap function from game v.
%
% Usage: [g bv lv]=Gap(v)
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
%   08/10/2010        0.1 beta        hme
%   09/14/2012        0.2             hme
%   10/27/2012        0.3             hme
%                



N=length(v);
[~, n]=log2(N);
% upper vector
Si=zeros(n,1);
k=1:n;
Si=bitset(N,k,0);
bv=v(N)-v(Si);


% Computing the gap function w.r.t. v.
g=zeros(1,N); % the gap vector w.r.t. v.
S=1:N;
Bm=bv(1); for ii=2:n, Bm=[Bm bv(ii) Bm+bv(ii)]; end
g=Bm-v;

a=cell(n,1);
lv=zeros(1,n);
% concession vector
for i=1:n
   a{i}=bitget(S,i)==1;
   lv(i)=min(g(a{i}));
end



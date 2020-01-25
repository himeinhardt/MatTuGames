function [ck kcvq]=k_convexQ(v)
% K_CONVEXQ checks whether the game v is k-convex.
%
%
% Usage: [ck kcvq]=k_convexQ(v)
% Define variables:
%  output:
%  ck       -- Returns 0 in case of a non k-convex game, otherwise 
%              a vector/list of numbers indicating which kind of
%              k-convexity has been discovered. 
%  kcvq     -- Returns a logical vector.
%  input:
%  v        -- A TU-game of length 2^n-1.
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



N=length(v);
[~, n]=log2(N);
vk_g=cell(n,1);
kcQ_g=zeros(n,1);

for k=1:n;
[vk_g{k} kcQ_g(k)]=k_cover(v,k);
gcQ(k)=convex_gameQ(vk_g{k});
end

kcvq=gcQ & kcQ_g'; 
J=1:n;
ck=J(kcvq);

if isempty(ck)
  ck=0;
else
end

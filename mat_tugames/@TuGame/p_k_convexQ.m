function [ck kcvq]=p_k_convexQ(clv)
% P_K_CONVEXQ checks whether the game v is k-convex.
% Using Matlab's PCT.
%
%
% Usage: [ck kcvq]=p_k_convexQ(clv)
%
% Define variables:
%  output:
%  ck       -- Returns 0 in case of a non k-convex game, otherwise 
%              a vector/list of numbers indicating which kind of
%              k-convexity has been discovered. 
%  kcvq     -- Returns a logical vector.
%
%  input:
%  clv        -- TuGame class object.
%

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

vk_g=cell(n,1);
kcQ_g=zeros(n,1);

parfor k=1:n;
[vk_g{k}, kcQ_g(k)]=k_cover(v,k);
gcQ(k)=convex_gameQ(vk_g{k});
end

kcvq=gcQ & kcQ_g'; 
J=1:n;
ck=J(kcvq);

if isempty(ck)
  ck=0;
else
end

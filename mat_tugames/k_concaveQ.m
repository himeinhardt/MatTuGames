function [ck kcvq]=k_concaveQ(v,tol)
% K_CONCAVEQ checks whether the game v is k-concave.
%
%
% Usage: [ck kcvq]=k_concaveQ(v)
% Define variables:
%  output:
%  ck       -- Returns 0 in case of a non k-concave game, otherwise 
%              a vector/list of numbers indicating which kind of
%              k-concaveity has been discovered. 
%  kcvq     -- Returns a logical vector.
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- A tolerance value. Default is set to tol=2*10^6*eps;
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/18/2016        0.8             hme
%                

if nargin<2
   tol=2*10^6*eps;
end


N=length(v);
[~, n]=log2(N);
vk_g=cell(n,1);
kcQ_g=zeros(n,1);

for k=1:n;
[vk_g{k} kcQ_g(k)]=k_anticover(v,k);
gcQ(k)=concave_gameQ(vk_g{k},tol);
end
kcvq=gcQ & kcQ_g'; 
J=1:n;
ck=J(kcvq);

if isempty(ck)
  ck=0;
else
end

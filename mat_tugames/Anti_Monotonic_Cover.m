function vc=Anti_Monotonic_Cover(v)
% ANTI_MONOTONIC_COVER computes the anti-monotonic cover of game v 
% or from excess game ex_v. 
% 
% Usage: v=Anti-Monotonic_Cover(v)
%

% Define variables:
%  output:
%  vc       -- Returns the anti-monotonic cover of game v 
%              or from excess game ex_v.
%
%  input:
%  v        -- A Tu-Game v or an excess game ex_v of length 2^n-1. 
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/28/2021        1.9             hme
%                

N=length(v);
[~, n]=log2(N);
vc=zeros(1,N);
for S=1:N;
   sS=SuperSets(S,n);
   vc(S)=min(v(sS));
end



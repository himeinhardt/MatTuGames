function vc=monotonic_cover(v);
% MONOTONIC_COVER computes the monotonic cover of game v 
% or from excess game ex_v. 
% 
% Usage: v=monotonic_cover(v)
%

% Define variables:
%  output:
%  vc       -- Returns the monotonic cover of game v 
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
%   06/20/2013        0.4             hme
%                

N=length(v);
[~, n]=log2(N);
vc=zeros(1,N);
for S=1:N;
   sS=SubSets(S,n);
   vc(S)=max(v(sS));
end



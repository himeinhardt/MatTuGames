function decQ=DecomposableQ(v,cs,tol)
% DECOMPOSABLEQ checks whether the game v is decomposable w.r.t. the coalition structure cs.
%
% Usage: DecQ=DecomposableQ(v,cs)
%
% Define variables:
%  decQ     -- Returns 1 (true) whenever the game v is decomposable w.r.t. cs, otherwise 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  tol      -- Tolerance value. Its default value is set to 10^7*eps.
%
% Example:
%  Define a game by
%  v=[2    5    5    7    8    9    9    2    4    7    7    9   10   11   11    3    5    8    8   10   11   12   12    4    6    9    9   11   12   13   13];
%     
% Choose a coalition structure cs=[7 24], then invoke    
%
%  decQ=DecomposableQ(v,cs)    
%   
%  that returns   
%  decQ_v = 1  
%
    
%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/25/2021        1.9             hme
%

if nargin<3    
   tol=10^7*eps;
end    
    
N=length(v);
lcs=length(cs);
dec=false(1,N);
sS=1:N;
for S=1:N    
    dec(S)=abs(sum(v(sS(ismember(sS,bitand(cs,S)))))-v(S))<tol;
end
decQ=all(dec);

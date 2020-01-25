function subg=SubGame(v,S)
% SUBGAME computes from (v,S) a subgame vS on S of the TU-game v.
%
% Usage: subg=SubGame(v,S)
% Define variables:
%  output:
%  subg     -- The corresponding sub-game of v on S.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  S        -- A coalition/set identified by its unique integer representation.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/26/2010        0.1 beta        hme
%   10/27/2012        0.3             hme
%                

N=length(v);
[~, n]=log2(N);
subS=SubSets(S,n);
subg=v(subS);

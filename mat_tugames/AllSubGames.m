function subGs=AllSubGames(v)
% ALLSUBGAMES computes all subgames of the TU-game v.
% 
%  
% Usage: subGs=AllSubGames(v)
% Define variables:
%  output:
%  subGs    -- The corresponding sub-game of v on S.
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
%                


N=length(v);
subGs=cell(1,N);
for k=1:N
  subGs{1,k}=SubGame(v,k);
end

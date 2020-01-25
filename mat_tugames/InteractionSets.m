function ias=InteractionSets(S,hs,n)
% INTERATCTIONSETS computes the interaction set for the position value. 
%
% Usage: ias=InteractionSets(S,hs,n)
%
% Define variables:
%  output:
%  ias      -- a system of interaction sets.
%
%  input:
%  S        -- A coalition S.
%  hs       -- A  hypergraph communication system like [1 7 14 15] 
%              for {[1],[1 2 3], [2 3 4],[1 2 3 4]}.
%  n        -- numbers of player involved in the game.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/30/2013        0.4             hme
%   05/15/2014        0.5             hme
% 

J=1:n;    
uc=2.^(J-1);    
sb=SubSets(S,n);
hs=sort(hs);
iuc=sb(ismembc(sb,uc));
slc=sb(ismembc(sb,hs));
ias=sort([iuc slc]);

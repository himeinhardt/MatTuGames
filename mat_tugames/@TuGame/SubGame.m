function subg=SubGame(clv,S)
% SUBGAME computes from (clv,S) a subgame vS on S of the TU-game v.
%
% Usage: subg=SubGame(clv,S)
%
% Define variables:
%  output:
%  subg     -- The corresponding sub-game of v on S.
%
%  input:
%  clv      -- TuGame class object.
%  S        -- A coalition/set identified by its unique integer representation.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/28/2012        0.3             hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
subS=SubSets(S,n);
subg=v(subS);

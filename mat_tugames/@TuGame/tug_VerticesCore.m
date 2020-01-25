function SOL=tug_VerticesCore(clv)
% TUG verifies game solution with the Mathematica Package TuGames.
%
% Define variables:
%  output:
%  core_vert  -- Matrix of core vertices. Output is numeric or a string.
%  crst       -- The core constraints.
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
%   11/02/2012        0.3             hme
%

n=clv.tuplayers;

math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Game to Mathematica ...')
w=gameToMama(clv);
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
disp('Computing the Vertices ...')
math('vert=VerticesCore[ExpGame]');
V=math('vert');
% unbalanced list
vert_v=math('math2matlab','vert[[1]]');
SOL=struct('Vertices',vert_v,'MVertices',V);
math('quit')

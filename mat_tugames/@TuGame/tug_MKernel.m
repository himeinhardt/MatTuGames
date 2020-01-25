function SOL=tug_MKernel(clv)
% TUG_MKernel computes a kernel element with the Mathematica Package TuGames.
%
% Usage: SOL=tug_MKernel(v)
% Define variables:
%  output:
%  SOL        -- Returns a kernel element based on Owen's algorithm (1974).
%                Field variable gives result in Matlab and Mathematica format.  
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
ker=math('mker=ModifiedKernel[ExpGame]');
ker_v=math('math2matlab','mker');
SOL=struct('Kernel',ker_v,'MKernel',ker);
math('quit')

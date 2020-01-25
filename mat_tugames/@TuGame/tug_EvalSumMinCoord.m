function SOL=tug_EvalSumMinCoord(clv)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_EvalSumMinCoord(v)
% Define variables:
%  output:
%  SOL        -- Returns the minimum sum of at most (n-1) inequalities of the unanimity coordinates constraints having nonnegative sums.
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
%   11/01/2012        0.3             hme
%

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;


math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Cost Game to Mathematica ...')
w=gameToMama(clv);
math('matlab2math','n1',n);
math('matlab2math','mw',w);
math('mv=Flatten[mw,1]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=PrependTo[mv,0];}');
math('ExpGame:=(DefineGame[T,mg];);');
mhd=math('muc=EvalSumMinCoord[ExpGame]');
svg=math('math2matlab','muc');
SOL=struct('EvalSumMinCoord',svg,'MEvalSumMinCoord',mhd);
math('quit')

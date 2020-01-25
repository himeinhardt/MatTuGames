function SOL=tug_ConvexUnanConditionQ(clv)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_ConvexUnanConditionQ(v)
% Define variables:
%  output:
%  SOL        -- Checking the sufficient Unanimity Conditions for Convexity
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/04/2012        0.3             hme
%

v=clv.tuvalues;
n=clv.tuplayers;

math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Game to Mathematica ...')
w=gameToMama(v);
math('matlab2math','n1',n);
math('matlab2math','mw',w);
math('mv=Flatten[mw,1]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=PrependTo[mv,0];}');
math('ExpGame:=(DefineGame[T,mg];);');
disp('Checking the sufficient Unanimity Conditions for Convexity ...');
msvg=math('mbv=ConvexUnanConditionQ[ExpGame]');
SOL=struct('ConvexUnanimityCondQ',msvg);
math('quit')

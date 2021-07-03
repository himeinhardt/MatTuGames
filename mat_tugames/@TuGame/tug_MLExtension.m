function SOL=tug_MLExtension(clv)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=clv.tug_ShapleyValueML
% Define variables:
%  output:
%  SOL        -- Returns the multi-linear extension of game v.  
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/26/2014        0.5             hme
%   07/02/2021        1.9             hme
%

n=clv.tuplayers;

math('quit')
pause(1)
math('$Version')
try 
    math('{Needs["TUG`"] }'); 
catch 
    math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }'); 
end
disp('Passing Game to Mathematica ...')
w=clv.gameToMama;
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
F=math('F=MLExtension[ExpGame]');
%f=math('math2matlab','F');
SOL=struct('MLEx','','Mmlex',F);
math('quit')

function SOL=tug_LowerSetQ(clv)
% TUG_LOWERSETQ verifies if the lower set exists. It calls the Mathematica Package TuGames.
% This function requires the Mathematica Symbolic Toolbox. It is available under the URL: 
% http://www.mathworks.com/matlabcentral/fileexchange/6044-mathematica-symbolic-toolbox-for-matlab-version-2-0
%   
% Usage: SOL=tug_LowerSetQ(v)   
% Define variables:
%  SOL          -- Returns 'True' or 'False' in Mathematica output format.
%                  It is a field variable.
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
w=gameToMama(clv);
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
crQ=math('LowerSetQ[ExpGame]');
SOL=struct('LowerSetQ',crQ);
math('quit')

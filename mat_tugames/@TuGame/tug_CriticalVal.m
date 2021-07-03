function SOL=tug_CriticalVal(clv)
% TUG_CriticalVal verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_CriticalVal(v)
% Define variables:
%  output:
%  SOL          -- A critical value to obtain a strong epsilon core that contains the kernel.
%                  Field variable gives result in Matlab and Mathematica format.
%
%  input:
%  clv          -- TuGame class object.
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

v=clv.tuvalues;
N=clv.tusize;
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
math('stx=Flatten[x1]');
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
fcv=math('mfcv=FirstCriticalVal[ExpGame]');
scv=math('mscv=SecondCriticalVal[ExpGame]');
tcv=math('mtcv=ThirdCriticalVal[ExpGame]');
fthcv=math('mfthcv=FourthCriticalVal[ExpGame]');
stcv=math('mstcv=StarCriticalVal[ExpGame]');
sstcv=math('msstcv=SecondStarCriticalVal[ExpGame]');
tstcv=math('mtstcv=ThirdStarCriticalVal[ExpGame]');
ffcv=math('mffcv=FifthCriticalVal[ExpGame]');
SOL=struct('FirstCriticalVal',fcv,'SecondCriticalVal',scv,'ThirdCriticalVal',tcv,'FourthCriticalVal',fthcv,'FifthCriticalVal',ffcv,'StarCriticalVal',stcv,'SecondStarCriticalVal',sstcv,'ThirdStarCriticalVal',tstcv);
math('quit')

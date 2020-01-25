function SOL=tug_CriticalVal(v)
% TUG_CriticalVal verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_CriticalVal(v)
% Define variables:
%  output:
%  SOL          -- A critical value to obtain a strong epsilon core that contains the kernel.
%                  Field variable gives result in Matlab and Mathematica format.
%  input:
%  v            -- A Tu-Game v of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/06/2011        0.1 beta        hme
%

% Here we assume that the user has represented the game correctly.
if nargin<1
    error('At least the game must be given!');
elseif nargin<2
N=length(v);
gr=dec2bin(N);
n=length(gr);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
else
    N=length(v);
    gr=dec2bin(N);
    n=length(gr);
    if (2^n-1)~=N
       error('Game has not the correct size!');
    end
end



math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Game to Mathematica ...')
w=gameToMama(v);
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

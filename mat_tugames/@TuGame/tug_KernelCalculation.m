function SOL=tug_KernelCalculation(clv)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_KernelCalculation(v)
% Define variables:
%  output:
%  SOL        -- A kernel element of game v and a list of kernel candidates.
%                Field variable gives result in Matlab and Mathematica format.
%
%  input:
%  clv      -- TuGame class object.
%
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
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
disp('Computing the Kernel ...')
math('sceps=SecondCriticalVal[ExpGame]');
math('{solker,pay}=KernelCalculation[ExpGame,EpsilonValue->sceps[[1,2]]]');
ker_el=math('solker');
cand=math('pay');
ker_cand=math('math2matlab','pay');
ker01_v=math('math2matlab','solker');
SOL=struct('Kernel',ker01_v,'MKernel',ker_el,'KernelCandidates',ker_cand,'MKernelCandidates',cand);
math('quit')

function SOL=tug_FindPreKernel(clv,y)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_FindPreKernel(v,y)
% Define variables:
%  output:
%  SOL        -- Returns a pre-kernel element of length (1xn).
%                For details see U. Faigle, W. Kern and J. Kuipers (1998).
%                Field variable gives result in Matlab and Mathematica format.
%
%  input:
%  clv        -- TuGame class object.
%  y          -- A payoff vector/matrix of length (1xn).
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
if nargin < 2
   if isa(clv,'TuSol')
      y=clv.tustpt;
   elseif isa(clv,'p_TuSol')
      y=clv.tustpt;
   else
      y=(v(N)/n)*ones(1,n);
   end
   if isempty(y)
     y=(v(N)/n)*ones(1,n);
   end
end

math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Game to Mathematica ...')
w=gameToMama(clv);
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('matlab2math','x1',y);
math('bds=Flatten[n1][[1]]');
math('stx=Flatten[x1]');
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
math('rtx=Rationalize[stx]');
disp('Computing the Pre-Kernel ...')
pk1=math('pk01=FindPreKernelSolution[ExpGame,rtx]');
prk_v=math('math2matlab','pk01');
mprkQ=math('prkQ=PreKernelQ[ExpGame,pk01]');
SOL=struct('PreKernel',prk_v,'MPreKernel',pk1, 'MPreKernel_Q',mprkQ);
math('quit')

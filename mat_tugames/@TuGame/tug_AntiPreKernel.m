function SOL=tug_AntiPreKernel(clv,y)
% TUG verifies game solution with the Mathematica Package TuGames.
%
% Define variables:
%  output:
%  core_vert  -- Matrix of core vertices. Output is numeric or a string.
%  crst       -- The core constraints.
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
      y=clv.tu_aprk;
   elseif isa(clv,'p_TuSol')
      y=clv.tu_aprk;
   else
      y=(v(N)/n)*ones(1,n);
   end
   if isempty(y)
     y=(v(N)/n)*ones(1,n);
   end
end

ly=length(y);
if n~=ly
  error('Vector has not the correct dimension!');
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
disp('Computing the Anti Pre-Kernel ...')
pk1=math('pk01=AntiPreKernelSolution[ExpGame,rtx,SolutionExact->False]');
prk_v=math('math2matlab','pk01');
SOL=struct('AntiPreKernel',prk_v,'MAntiPreKernel',pk1);
math('quit')

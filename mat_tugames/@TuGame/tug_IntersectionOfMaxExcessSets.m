function SOL=tug_IntersectionOfMaxExcessSets(clv,y)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_IntersectionOfMaxExcessSets(v,y)
% Define variables:
%  output:
%  SOL        -- Determines if the set of proper coalitions with maximum surpluses has an empty intersection.
%                It is a field variable.
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
math('vec=Table[0,{i,bds}]');
math('vec[[1]]=v[T]');
math('vec');
math('perm=Permutations[vec]');
math('rtx=Rationalize[stx]');
mb1=math('mb01=IntersectionOfMaxExcessSets[ExpGame,rtx]');
SOL=struct('IntersectionOfMaxExcessSets',mb1);
math('quit')

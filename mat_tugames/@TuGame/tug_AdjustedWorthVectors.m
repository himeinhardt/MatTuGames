function SOL=tug_AdjustedWorthVectors(clv,int_k)
% TUG_LargestAmount verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_AdjustedWorth(v,int_k)
% Define variables:
%  output:
%  SOL        -- Adjusted marginal worth vector. 
%                Field variable gives result in Matlab and Mathematica format.
%
%  input:
%  clv        -- TuGame class object.
%  int_k      -- An integer samller than or equal to n.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/01/2012        0.3             hme
%

n=clv.tuplayers;

if nargin<2
  error('An integer not larger than n must be given!')
else
   if isnumeric(int_k)
      if length(int_k)>1
        error('Second argument must be an integer not larger than n!')
      end
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
math('matlab2math','intk',int_k);
math('k=Rationalize[Flatten[intk,1][[1]]]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('expg=Flatten[mg1,1]');
math('{T,mg=PrependTo[expg,0]}');
math('ExpGame:=(DefineGame[T,mg];);');
madw=math('adw=AdjustedWorthVectors[ExpGame,k]');
adw_v=math('math2matlab','adw');
SOL=struct('AdjustedWorth',adw_v,'MAdjustedWorth',madw);
math('quit')

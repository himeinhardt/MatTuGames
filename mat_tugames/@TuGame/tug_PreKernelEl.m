function SOL=tug_PreKernelEl(clv,y)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_PreKernelEl(v,y)
% Define variables:
%  output:
%  SOL        -- Returns a pre-kernel point based on Algorithm 6.3 by Meinhardt 2010.
%                Field variable gives result in Matlab and Mathematica format.
%
%  input:
%  clv        -- TuGame class object.
%  y          -- A payoff vector/matrix of length (1xn) or size (mxn).
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
try 
    math('{Needs["TUG`"] }'); 
catch 
    math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }'); 
end
disp('Passing Game to Mathematica ...')
w=gameToMama(clv);
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('matlab2math','x1',y);
math('bds=Flatten[n1][[1]]');
if length(y)>1
 if isvector(y)
    math('matlab2math','x1',y);
    math('stx=Flatten[x1]');
 elseif ismatrix(y)
    math('matlab2math','x1',y);
    math('stx=x1');
 else
 end
else
  disp('Payoff x is neither a vector nor a matrix!')
  disp('Falling back on default vector!')
  y=(v(N)/n)*ones(1,n);
  math('matlab2math','x1',y);
  math('stx=Flatten[x1]');
end
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
math('rtx=Rationalize[stx]');
disp('Computing the Pre-Kernel ...')
pk1=math('pk01=PreKernelElement[ExpGame,rtx, SmallestCardinality -> True]');
prk_v=math('math2matlab','pk01');
mprkQ=math('prkQ=PreKernelQ[ExpGame,pk01]');
SOL=struct('PreKernel',prk_v,'MPreKernel',pk1, 'MPreKernel_Q',mprkQ);
math('quit')

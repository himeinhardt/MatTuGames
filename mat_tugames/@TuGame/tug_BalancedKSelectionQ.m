function SOL=tug_BalancedKSelectionQ(clv,y,k)
% TUG_BALANCEDKSELECTION verifies if the induced collections are K-balanced.
%
% Usage: SOL=tug_BalancedKSelectionQ(v,y,k)
% Define variables:
%  output:
%  SOL        -- Returns 'True' or 'False' in Mathematica output format.
%                It is a field variable. 
%
%  input:
%  clv        -- TuGame class object.
%  y          -- A payoff vector/matrix of length (1xn) or size (mxn).
%  k          -- An integer between 2 and n.
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
k=3;
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
elseif nargin < 3
  k=3;
else
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
w=gameToMama(v);
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('matlab2math','k1',k);
math('bds=Flatten[n1][[1]]');
math('k=Rationalize[Flatten[k1][[1]]]');
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
py=math('rtx=Rationalize[stx]');
disp('Is the selection K-balanced?...')
mprkQ=math('bsQ=SelectionKBalancedQ[ExpGame,rtx,k,Silent->True,Tight->True]');
SOL=struct('BalancedKSelectionQ',mprkQ);
math('quit')

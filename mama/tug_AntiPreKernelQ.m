function SOL=tug_AntiPreKernelQ(v,y)
% TUG verifies game solution with the Mathematica Package TuGames.
%
% Define variables:
%  output:
%  core_vert  -- Matrix of core vertices. Output is numeric or a string.
%  crst       -- The core constraints.
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%  y          -- A payoff vector/matrix of length (1xn).
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
    y=(v(N)/n)*ones(1,n);
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
disp('Is payoff a Pre-Kernel element?...')
mprkQ=math('prkQ=AntiPreKernelQ[ExpGame,rtx]');
SOL=struct('Payoff',y,'MPayoff',py, 'MAntiPreKernel_Q',mprkQ);
math('quit')

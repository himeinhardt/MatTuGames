function SOL=tug_EqClass(v,y)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_EqClass(v,y)
% Define variables:
%  output:
%  SOL        -- Computes the equivalence class from the set of most effective coalitions.
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%  y          -- A payoff vector of length (1xn).
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
    smc=1;
else
    N=length(v);
    gr=dec2bin(N);
    n=length(gr);
    ly=length(y);
    if (2^n-1)~=N
       error('Game has not the correct size!');
    end
    if n~=ly
       error('Vector has not the correct dimension!');
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
math('matlab2math','x1',y);
math('bds=Flatten[n1][[1]]');
math('stx=Flatten[x1]');
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
math('rtx=Rationalize[stx]');
stv1=math('eqc=ImputationToEqClass[ExpGame,rtx,EffVector->False]');
SOL=struct('EqClass',stv1);
math('quit')

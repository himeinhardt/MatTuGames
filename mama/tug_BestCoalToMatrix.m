function SOL=tug_BestCoalToMatrix(v,y)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_BestCoalToMatrix(v,y) 
% Define variables:
%  output:
%  SOL        -- Returns a set of equivalence class as a matrix.
%                Field variable gives result in Matlab and Mathematica format.
%
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%  y            -- A payoff vector of length (1xn).
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
disp('Computing the Pre-Kernel ...')
bsc1=math('bsc01=BestCoalitions[ExpGame,rtx]');
math('{Q,Qh,subE}=BestCoalToMatrix[bsc01,T]');
mQ1=math('math2matlab','Q[[1]]');
mQ3=math('math2matlab','{Q[[3]]}');
SOL=struct('MatrixQ',mQ1,'DetQ',mQ3,'BestCoalitions',bsc1);
math('quit')

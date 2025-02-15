function SOL=tug_EvalSumMinCoord(v)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_EvalSumMinCoord(v)
% Define variables:
%  output:
%  SOL        -- Returns the minimum sum of at most (n-1) inequalities of the unanimity coordinates constraints having nonnegative sums.
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/06/2011        0.1 beta        hme
%   07/02/2021        1.9             hme
%

% Here we assume that the user has represented the game correctly.

if nargin<1
    error('At least the cost game must be given!');
elseif nargin<2
N=length(v);
gr=dec2bin(N);
n=length(gr);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
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
try 
    math('{Needs["TUG`"] }'); 
catch 
    math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }'); 
end
disp('Passing Cost Game to Mathematica ...')
w=gameToMama(v);
math('matlab2math','n1',n);
math('matlab2math','mw',w);
math('mv=Flatten[mw,1]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=PrependTo[mv,0];}');
math('ExpGame:=(DefineGame[T,mg];);');
mhd=math('muc=EvalSumMinCoord[ExpGame]');
svg=math('math2matlab','muc');
SOL=struct('EvalSumMinCoord',svg,'MEvalSumMinCoord',mhd);
math('quit')

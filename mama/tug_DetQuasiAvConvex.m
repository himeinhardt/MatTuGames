function SOL=tug_DetQuasiAvConvex(v)
% TUG verifies game solution with the Mathematica Package TuGames.
%
% Usage: SOL=tug_DetQuasiAvConvex(v)
% Define variables:
%  output:
%  SOL        -- Returns a quasi average game from game v.
%                It is a field variable.
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
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Game to Mathematica ...')
w=gameToMama(v);
math('matlab2math','n1',n);
math('matlab2math','mw',w);
math('mv=Flatten[mw,1]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=PrependTo[mv,0];}');
math('ExpGame:=(DefineGame[T,mg];);');
disp('Determing an Average Convex Game ...');
msvg=math('mbv=DetQuasiAvConvex[ExpGame]');
svg=math('math2matlab','mbv');
svg(:,1)=[];
sv_g=gameToMatlab(svg);
SOL=struct('TuGame',sv_g,'MTuGame',msvg);
math('quit')

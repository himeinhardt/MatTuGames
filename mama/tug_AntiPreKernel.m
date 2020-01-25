function SOL=tug_AntiPreKernel(v,y)
% TUG verifies game solution with the Mathematica Package TuGames.
%
% Define variables:
%  output:
%  core_vert  -- Matrix of core vertices. Output is numeric or a string.
%  crst       -- The core constraints.
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%  y            -- A payoff vector/matrix of length (1xn).
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
disp('Computing the Anti Pre-Kernel ...')
pk1=math('pk01=AntiPreKernelSolution[ExpGame,rtx]');
prk_v=math('math2matlab','pk01');
SOL=struct('AntiPreKernel',prk_v,'MAntiPreKernel',pk1);
math('quit')

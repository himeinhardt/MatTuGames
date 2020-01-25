function SOL=tug_CddGmpVerticesCore(v)
% TUG verifies game solution with the Mathematica Package TuGames.
%
% Define variables:
%  output:
%  core_vert  -- Matrix of core vertices. Output is numeric or a string.
%  crst       -- The core constraints.
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%  method     -- A string to call a method from the cdd-library.
%                Permissible methods are:
%                'float' that is, results are given by real numbers.
%                Default is 'float'.
%                'gmp' that is, results are given by rational numbers.
%                Choose this method whenever the result with 'float'
%                is not as expected. This method needs more time to complete.
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
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('Print["mg=",mg];');
math('ExpGame:=(DefineGame[T,mg];);');
disp('Computing the Vertices ...')
math('vert=CddGmpVerticesCore[ExpGame]');
V=math('vert');
vert_v=math('math2matlab','vert');
SOL=struct('Vertices',vert_v,'MVertices',V);
math('quit')

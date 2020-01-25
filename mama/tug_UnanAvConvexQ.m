function SOL=tug_UnanAvConvexQ(hd)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_UnanAvConvexQ(hd)
% Define variables:
%  output:
%  SOL         -- Checks if the coordinates satisfy the sufficient and necessary condition of average convexity
%                 in terms of unanimity coordinates.
%                 Returns 'True' or 'False'.
%                 Field variable gives result in Matlab and Mathematica format.
%  input:
%  hd          -- The unanimity coordinates of length 2^n-1.
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
    error('At least the unanimity coordinates must be given!');
elseif nargin<2
N=length(hd);
gr=dec2bin(N);
n=length(gr);
    if (2^n-1)~=N
      error('Unanimity coordninates have not the correct size!');
    end
else
    N=length(hd);
    gr=dec2bin(N);
    n=length(gr);
    if (2^n-1)~=N
       error('Unanimity coordninates have not the correct size!');
    end
end



math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing unanimity coordninates to Mathematica ...')
whd=gameToMama(hd);
math('matlab2math','n1',n);
math('matlab2math','mhd',whd);
math('uc=Flatten[mhd,1]');
math('bds=Flatten[n1][[1]]');
%math('Print["bds=",bds]');
math('T=Flatten[Range[n1]]');
math('{T,coord=PrependTo[uc,0];}');
%math('Print["T=",T]');
%math('Print["coord=",coord]');
msvg=math('mbv=UnanAvConvexQ[coord,T]');
SOL=struct('AverageConvexQ',msvg);
math('quit')

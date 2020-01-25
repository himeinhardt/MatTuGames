function SOL=tug_CharacteristicValues(hd)
% TUG verifies game/property solution with the Mathematica Package TuGames.
%
% Usage: SOL=tug_CharacteristicValues(hd)
% Define variables:
%  output:
%  SOL        -- A TU-game of length 2^n-1. 
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  hd         -- Unanimity coordinates hd of length 2^n-1.
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
math('T=Flatten[Range[n1]]');
math('{T,coord=PrependTo[uc,0];}');
msvg=math('mbv=CharacteristicValues[coord,T]');
w=math('math2matlab','mbv');
w(:,1)=[];
v=gameToMatlab(w);
SOL=struct('CharacteristicValues',v,'MCharacteristicValues',msvg);
math('quit')

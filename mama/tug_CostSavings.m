function SOL=tug_CostSaving(c)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_CostSaving(c)
% Define variables:
%  output:
%  SOL        -- A Tu-Game v of length 2^n-1.
%  input:
%  c          -- A cost game of length 2^n-1.
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
N=length(c);
gr=dec2bin(N);
n=length(gr);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
else
    N=length(c);
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
disp('Passing Cost Game to Mathematica ...')
math('matlab2math','n1',n);
math('matlab2math','c1',c);
math('cv=Flatten[c1,1]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mcv=PrependTo[cv,0];}');
disp('Determing Savings Game ...')
msvg=math('mbv=CostSaving[mcv,T]');
svg=math('math2matlab','mbv');
svg(:,1)=[];
sv_g=gameToMatlab(svg);
SOL=struct('SavingsGame',sv_g,'MSavingsGame',msvg);
math('quit')

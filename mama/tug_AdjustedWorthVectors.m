function SOL=tug_AdjustedWorth(v,int_k)
% TUG_LargestAmount verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_AdjustedWorth(v,int_k)
% Define variables:
%  output:
%  SOL        -- Adjusted marginal worth vector. 
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  v          -- A Tu-Game v of length 2^n-1.
%  int_k      -- An integer samller than or equal to n.

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
    if length(int_k)>1
      error('Second argument must be an integer not larger than n!')
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
math('matlab2math','intk',int_k);
math('k=Rationalize[Flatten[intk,1][[1]]]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('expg=Flatten[mg1,1]');
math('{T,mg=PrependTo[expg,0]}');
math('ExpGame:=(DefineGame[T,mg];);');
madw=math('adw=AdjustedWorthVectors[ExpGame,k]');
adw_v=math('math2matlab','adw');
SOL=struct('AdjustedWorth',adw_v,'MAdjustedWorth',madw);
math('quit')

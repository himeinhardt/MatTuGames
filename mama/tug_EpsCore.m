function SOL=tug_EpsCore(v)
% TUG_EPSCORE computes the epsilon value of the least core. It calls the Mathematica Package TuGames.
% This function requires the Mathematica Symbolic Toolbox. It is available under the URL: 
% http://www.mathworks.com/matlabcentral/fileexchange/6044-mathematica-symbolic-toolbox-for-matlab-version-2-0
% 
% Usage: SOL=tug_EpsCore(v)   
% Define variables:
%  output:
%  SOL          -- Returns the epsilon value of the least core and a payoff vector belonging to it.
%                  Field variable gives result in Matlab and Mathematica format.
%  input:
%  v            -- A Tu-Game v of length 2^n-1.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/08/2011        0.1 beta        hme
%   07/02/2021        1.9             hme
%

    
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
try 
    math('{Needs["TUG`"] }'); 
catch 
    math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }'); 
end
disp('Passing Game to Mathematica ...')
w=gameToMama(v);
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('{T,mg=FlattenAt[PrependTo[mg1,0],2];}');
math('ExpGame:=(DefineGame[T,mg];);');
crQ=math('EpsCore[ExpGame]');
SOL=struct('CorePayoff',crQ);
math('quit')

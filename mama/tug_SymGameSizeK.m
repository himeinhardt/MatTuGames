function SOL=tug_SymGameSizeK(n,k,val)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_SymGameSizeK(n,k,val)
%  Example: tug_SymGameSizeK(5,3,2)
% Define variables:
%  output:
%  SOL        -- Returns a symmetric game.
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  n          -- Persons involved to specify the grand coalition, must be an integer value greater than 3.
%  k          -- Specifies a cylce of coalitions of size k which should be assigned with the worth 'val'. 
%                Must be an integer value.
%  val        -- Value assinged to the coalitions. Integer value. Default is 3.
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
    error('The number of players for the player set must be specified! It should have at least 4-persons!');
    k=2;
    val=2;
elseif nargin<2
    error('The size of the subcoalition must be given by an integer smaller than n!');
    val=2;
else
end

%val
N=2^n-1;

if k>n
   error('The sub-coalition has size lager than n!');
end


math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Data to Mathematica ...')
math('matlab2math','t',n);
math('t1=Rationalize[Flatten[t,1][[1]]]');
math('T=Range[t1]');
math('matlab2math','mv1',val);
math('mv=Rationalize[Flatten[mv1,1][[1]]]');
math('matlab2math','k1',k);
math('k=Rationalize[Flatten[k1,1][[1]]]');
disp('Determing Symmetric Game with size k ...')
mabv=math('gsym=SymGameSizeK[T,k,mv]');
mbg=math('math2matlab','gsym');
mbg(:,1)=[];
sym_v=gameToMatlab(mbg);
SOL=struct('SymGame',sym_v,'MSymGame',mabv,'MGame',mbg);
math('quit')

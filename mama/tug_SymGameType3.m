function SOL=tug_SymGameType3(n,S,val)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_SymGameType3(n,S,val)
%  Example: SOL=tug_SymGameType3(5,22,4)
% Define variables:
%  output:
%  SOL        -- Returns a symmetric game.
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  n          -- Persons involved to specify the grand coalition, must be an integer value greater than 3.
%  S          -- Coalition of size not less than 3. Must be an integer value.
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
    S=7;
    val=3;
elseif nargin<2
    error('A subcoalition of size not less than 3 must be given!');
    val=3;
else
%    val=3;
end

N=2^n-1;
sS=dec2bin(S,n)-'0';

if S>=N
   error('The grand coalition must be larger than a sub-coalition!');
elseif nnz(sS)<3
   error('The sub-coalition has size smaller than 3!');
else
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
math('matlab2math','S1',S);
math('S1=Rationalize[Flatten[S1,1][[1]]]');
math('dS=Reverse[IntegerDigits[S1,2,t1]]');
math('S2=DeleteCases[dS*T,0]');
disp('Determing Symmetric Game Type 3 ...')
mabv=math('gsym=SymGameType3[T,S2,mv]');
mbg=math('math2matlab','gsym');
mbg(:,1)=[];
sym_v=gameToMatlab(mbg);
SOL=struct('SymGame',sym_v,'MSymGame',mabv,'MGame',mbg);
math('quit')

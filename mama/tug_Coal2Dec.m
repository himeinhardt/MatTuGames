function COAL=tug_Coal2Dec(n)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: COAL=tug_Coal2Dec(n)
% Define variables:
%  output:
%  COAL       -- The list of proper coalitions in Mathematica order.
%                Field variable gives result in Matlab and Mathematica format.
%
%  input:
%  n          -- The number of persons involved.
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
if nargin <1
n=3;
elseif length(n)>1
  error('Input must be an integer equal to or larger than 2!')
else 
end

math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Input to Mathematica ...')
math('matlab2math','n1',n);
math('n1=Flatten[n1,1][[1]]');
math('n2=Rationalize[n1]');
disp('Determing proper coalitions in decimal representation ...')
math('T=Range[n2]');
mcl=math('Subsets[T]');
mps=math('ps=Coal2Dec[n2]');
S=math('math2matlab','ps');
COAL=struct('Coalitions',S,'MProperCoalitions',mps,'MCoalitions',mcl);
math('quit')

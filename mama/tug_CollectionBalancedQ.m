function SOL=tug_CollectionBalancedQ(cS,n)
% TUG_COLLECTIONBALANCEDQ verifies if the collection of coalitions cS is balanced 
%
% Usage: SOL=tug_CollectionBalancedQ(cS,n)
% Define variables:
%  output:
%  SOL        -- Yields 'True' or 'False' in Mathematica output format.
%                It is a field variable.
%  input:
%  cS         -- A collection of coalitions payoff w.r.t. to n-persons.
%  n          -- An integer specifying the number of persons involved.

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



math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
math('matlab2math','t',n);
math('t1=Rationalize[Flatten[t,1][[1]]]');
math('T=Range[t1]');
math('matlab2math','S1',cS);
math('S=Rationalize[Flatten[S1,1]]');
math('cS=Reverse[IntegerDigits[#,2,t1]& /@ S]');
math('S2=DeleteCases[#*T,0] & /@ cS');
disp('Is the collection balanced?...')
mprkQ=math('bsQ=BalancedSelectionQ[S2,Tight->True]');
SOL=struct('CollectionBalancedQ',mprkQ);
math('quit')

function SOL=tug_TalmudicRule(E,d)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_TalmudicRule(E,d)
% Define variables:
%  output:
%  SOL        -- Returns a division rule of a generalized contested garment problem of size n.
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  E          -- Estate, it is an integer value.
%  d          -- Claimants vector of length (1xn).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/09/2011        0.1 beta        hme
%

% Here we assume that the user has represented the game correctly.
if nargin<1
    error('The Estate and the claims vector must be provided!');
elseif nargin<2
    error('The claims vector must be provided!');
else
  if isvector(d)==0
     error('The claims must be a vector of length greater or equal to 2!'); 
  else
  end
  if length(E)~=1
     error('The Estate must be an integer!');
  else
  end
end



math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Contested Garment Problem to Mathematica ...')
math('matlab2math','est',E);
math('est1=Flatten[est,1][[1]]');
math('matlab2math','d1',d);
math('clv=Flatten[d1,1]');
disp('Computing the Contested Garment Allocation...')
mcgr=math('cgr=TalmudicRule[est1,clv]');
cg_rl=math('math2matlab','cgr');
SOL=struct('TalmudicRule',cg_rl,'MTalmudicRule',mcgr);
math('quit')

function SOL=tug_WeightedMajority(th,wghs)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_WeightedMajority(th,wghs)
% Define variables:
%  output:
%  SOL        -- Returns a worth vector of a weighted majority game.
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  th         -- The quota, that is, an integer value.
%  wghs       -- A list of weights of length (1xn).
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
    error('The Estate and the claims vector must be provided!');
elseif nargin<2
    error('The claims vector must be provided!');
else
  if isvector(wghs)==0
     error('The weights must be a vector of length greater or equal to 2!');
  else
  end
  if length(th)~=1
     error('The threshold must be an integer!');
  else
  end
end
n=length(wghs);

math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing Bankruptcy Situation to Mathematica ...')
math('matlab2math','n1',n);
math('matlab2math','thv',th);
math('ml=Rationalize[Flatten[thv,1][[1]]]');
math('matlab2math','gw',wghs);
math('wgh=Flatten[gw,1]');
math('PrependTo[wgh,ml]')
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
disp('Determing Weighted Majority Game ...')
mabv=math('mbv=WeightedMajority[T,wgh]');
mbg=math('math2matlab','mbv');
mbg(:,1)=[];
bv=gameToMatlab(mbg);
SOL=struct('WeightedMajority',bv,'MWeightedMajority',mabv);
math('quit')

function SOL=tug_DetRandCoord(sed_val,n)
% TUG verifies game solution/property with the Mathematica Package TuGames.
%
% Usage: SOL=tug_DetRandCoord(sed_val,n)
% Define variables:
%  output:
%  SOL        -- Returns the unanimity coordinates.
%                Field variable gives result in Matlab and Mathematica format.
%  core_vert  -- Matrix of core vertices. Output is numeric or a string.
%  crst       -- The core constraints.
%
%  input:
%  sedval    -- An integer, for instance, 1234.
%  n          -- The number of persons involved.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/28/2013        0.5             hme
%

N=2^n-1;
bn=binomial(n,2);


math('quit')
pause(1)
math('$Version')
math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }');
disp('Passing input arguments to Mathematica ...')
math('matlab2math','n1',n);
math('matlab2math','sed',sed_val);
math('msed=ToExpression[Rationalize[First[sed[[1]]]]]');
math('matlab2math','mbn',bn);
math('SeedRandom[msed]');
math('coord=RandomInteger[20,mbn[[1]]];');
math('T=Flatten[Range[n1]]');
math('mbv=DetUCoord[coord,T]');
math('gvec=CharacteristicValues[mbv,T];');
svg=math('math2matlab','gvec');
svg(:,1)=[];
sv_g=gameToMatlab(svg);
SOL=struct('Game',sv_g,'MGame',svg);
math('quit')

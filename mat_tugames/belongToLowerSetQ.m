function blsQ=belongToLowerSetQ(v,x)
% BELONGTOLOWERSETQ checks if the imputation x belongs to the lower set.
% using MPT3.
%
%  Usage: busQ=belongToLowerSetQ(v,x) 
%
%
% Define variables:
%  output:
%  blsQ     -- Returns true (1) if x belongs to the lower set,
%              otherwise false (0).
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n)
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/20/2014        0.6             hme
%


lsQ=lowersetQ(v);
if lsQ==1
   vert=CddLowerSetVertices(v);
   x0=x';
   Pu = Polyhedron(vert);
   blsQ=Pu.contains( x0 );
else
  error('Lower set is empty!');
end

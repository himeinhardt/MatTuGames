function busQ=belongToUpperSetQ(v,x)
% BELONGTOUPPERSETQ  checks if the imputation x belongs to the upper set
% using MPT3.
%
%  Usage: busQ=belongToUpperSetQ(v,x) 
%
%
% Define variables:
%  output:
%  busQ     -- Returns true (1) if x belongs to the upper set,
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


usQ=uppersetQ(v);
if usQ==1
   vert=CddUpperSetVertices(v);
   x0=x';
   Pu = Polyhedron(vert);
   busQ=Pu.contains( x0 );
else
  error('Upper set is empty!');
end

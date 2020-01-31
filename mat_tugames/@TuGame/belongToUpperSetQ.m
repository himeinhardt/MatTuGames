function busQ=belongToUpperSetQ(clv,x)
% BELONGTOUPPERSETQ  checks if the imputation x belongs to the upper set
% using MPT3.
%
%  Usage: busQ=belongToUpperSetQ(clv,x) 
%
%
% Define variables:
%  output:
%  busQ     -- Returns true (1) if x belongs to the upper set,
%              otherwise false (0).
%  input:
%  clv      -- TuGame class object.
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


usQ=clv.uppersetQ;
if usQ==1
   vert=clv.CddUpperSetVertices;
   x0=x';
   Pu = Polyhedron(vert);
   busQ=Pu.contains( x0 );
else
  error('Upper set is empty!');
end

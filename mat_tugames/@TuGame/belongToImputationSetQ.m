function bisQ=belongToImputationSetQ(clv,x)
% BELONGTOIMPUTATIONSETQ checks if the payoff x belongs to the imputation set
% using MPT3.
%
%  Usage: bisQ=belongToImputationSetQ(clv,x) 
%
%
% Define variables:
%  output:
%  bisQ     -- Returns true (1) if x belongs to the imputation set,
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

n=clv.tuplayers;
N=clv.tusize;
v=clv.tuvalues;

esQ=clv.tuessQ;
if esQ==1
   vert=clv.CddImputationVertices;
   x0=x';
   Pu = Polyhedron(vert);
   bisQ=Pu.contains( x0 );
else
  error('Imputation set is empty!');
end

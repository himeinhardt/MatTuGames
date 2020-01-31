function bisQ=belongToImputationSetQ(v,x)
% BELONGTOIMPUTATIONSETQ checks if the payoff x belongs to the imputation set
% using MPT3.
%
%  Usage: bisQ=belongToImputationSetQ(v,x) 
%
%
% Define variables:
%  output:
%  bisQ     -- Returns true (1) if x belongs to the imputation set,
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



N=length(v);
[~, n]=log2(N);

k=1:n;
vi=v(bitset(0,k));
esQ=sum(vi)<=v(N);
if esQ==1
   vert=CddImputationVertices(v);
   x0=x';
   Pu = Polyhedron(vert);
   bisQ=Pu.contains( x0 );
else
  error('Imputation set is empty!');
end

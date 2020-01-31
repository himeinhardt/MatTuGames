function blcQ=CddBelongToLeastCoreQ(clv,x)
% CDDBELONGTOLEASTCOREQ checks if the imputation x belongs to the least core
% using MPT3.
%
%  Usage: blcQ=CddBelongToLeastCoreQ(clv,x)
%
%
% Define variables:
%  output:
%  blcQ     -- Returns true (1) if x belongs to the least core,
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


LC=clv.CddLeastCoreVertices();
x0=x';
Plc = Polyhedron(LC.crv);
blcQ=Plc.contains( x0 );

function blcQ=CddBelongToLeastCoreQ(v,x)
% CDDBELONGTOLEASTCOREQ checks if the imputation x belongs to the least core
% using MPT3.
%
%  Usage: blcQ=CddBelongToLeastCoreQ(v,x)
%
%
% Define variables:
%  output:
%  blcQ     -- Returns true (1) if x belongs to the least core,
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


LC=CddLeastCoreVertices(v);
x0=x';
Plc = Polyhedron(LC.crv);
blcQ=Plc.contains( x0 );

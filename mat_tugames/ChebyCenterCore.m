function ChC=ChebyCenterCore(v) 
% CHEBYCENTERCORE computes the the Cheby Center of the core of game v using MPT3.
%
%  Usage: ChC=ChebyCenterCore(v)
%
%
% Define variables:
%  Structure elements of BTSGQ
%  output:
%  ChepyC   -- The Chepy Center of the core of game v if the core exists,
%              otherwise an inf-vector of length n.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/25/2021        1.9.1             hme
%


N=length(v);
[~, n]=log2(N);
try 
   crQ=CddCoreQ(v);
catch
   crQ=coreQ(v);
end	

if crQ==0
   ChC=inf(1,n);
else
   crv=CddCoreVertices(v);
   Pspc = Polyhedron(crv);
   Pspc.computeHRep();
   chC=Pspc.chebyCenter();
   ChC=chC.x';
end	

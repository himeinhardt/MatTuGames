function wbsQ=belongToWeberSetQ(v,x)
% BELONGTOWEBERSETQ  checks if the imputation x belongs to the Weber set
% using MPT3.
%
%  Usage: busQ=belongToWeberSetQ(v,x) 
%
%
% Define variables:
%  output:
%  wbsQ     -- Returns true (1) if x belongs to the Weber set,
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
%   02/12/2022        1.9.1             hme
%


vert=CddWeberSet(v);
x0=x';
Pwb = Polyhedron(vert);
wbsQ=Pwb.contains( x0 );

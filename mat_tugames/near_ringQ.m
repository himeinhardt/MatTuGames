function nrQ=near_ringQ(sC,n)
% NEAR_RINGQ checks whether the collection cS is a near ring.
%
% Example: 
%     sC=[2    7   11   13   14];
%     nrQ=near_ringQ(sC,4), returns true, that is, 
%     1.
%
% Usage: nrQ=near_ringQ(sC,n)
% Define variables:
%  output:
%  nrQ       -- Returns true (1) is the collection sC is a near
%               ring, otherwise false (0).
%
%  input:
%  sC       -- An array of intergers, which represents a unique 
%              integer representation of the collection sC.
%  n        -- A positive number, which indicates the number of 
%              players in a Tu-game

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/30/2014        0.6             hme
%

lsc=length(sC);
N=2^n-1;
sC1=unique([sC,N]);

for k=1:lsc
    sor(k)=isempty(setdiff(unique(bitor(sC(k),sC1)),sC1));
    sad(k)=all(setdiff(unique(bitand(sC(k),sC1)),sC1)==0);
end
nrQ = all(sor | sad);
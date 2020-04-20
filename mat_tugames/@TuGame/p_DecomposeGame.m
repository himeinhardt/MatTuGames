function DecG=DecomposeGame(clv)
% P_DECOMPOSEGAME computes the unique decomposition of a TU-game v
% from a linear basis using Matlab's PCT. The direct sum of w and z is equal to v.
% For n>12 this function needs some time to complete.
%
% Usage: DecG=p_DecomposeGame(clv)
%
% Define variables:
%  output:     structure elements
%  w        -- A TU-game from vector space W.
%  z        -- A TU-game from vector space Z.
%              Notice that the game space G is the direct sum
%              of the spaces W and Z.
%  input:
%  clv      -- TuGame class object.
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/11/2014        0.5             hme
%

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

if n==1
   w=v;z=[];
   DecG=struct('w',w,'z',z);
   return; 
end    

[cfc,lb]=clv.p_coeff_linearbasis('sparse');
k=1:n;
sC=2.^(k-1);
Z=lb(:,sC);
cz=cfc(sC);
z=cz*Z';
tv=true(1,N);
tv(sC)=false;
W=lb(:,tv);
cw=cfc(tv);
w=cw*W';
DecG=struct('w',w,'z',z);

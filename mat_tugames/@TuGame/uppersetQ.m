function usQ=uppersetQ(clv)
% UPPERSETQ checks if the lower set is non empty.
% Is a kernel catcher under certain conditions if it is non-empty.
%
%  Usage: usQ=uppersetQ(v) 
%
%
% Define variables:
%  output:
%  usQ      -- Returns true (1) if the upper set is non empty,
%              otherwise false (0).
%  input:
%  clv      -- TuGame class object.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/05/2014        0.6             hme
%
v=clv.tuvalues;
N=clv.tusize;
pa=clv.proper_amount();
usQ=sum(pa)>=v(N);

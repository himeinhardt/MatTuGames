function w=gameToMama(v)
%  GAMETOMAMA converts a TU-game v that is based on 
%  a unique integer representation (like Matlab)
%  into a game w that is based on a generic power set 
%  representation (like Mathematica).
%
%  Example: 
%   A TU-game v that is based on the following binary order of coalitions
%   S=[1 2 12 3  13   23 123  4  14  24 124  34 134 234 1234] (full form)
%   whereas the same set of coalitions characterized in terms of their unique 
%   integer representation is given by
%   S=[1 2  3 4   5    6  7   8   9  10  11  12  13  14   15].
%   Now, the following coalitions values of game
%   v=[0 0 1  0  1/2  1  5/4  0  1   0   1   1   1   1   2],
%   will be converted to the TU-game w given by
%   w=[0 0 0 0 1  1/2  1  1  0  1  5/4  1   1   1   2 ],
%   based on the following power set representation as given by
%   S=[1 2 3 4 12 13  14 23 24 34 123 124 134 234 1234]. 
%   Use function clToMatlab() to verify these results.
%
% Usage: w=gameToMama(v)
% Define variables:
%  output:
%  w        -- A Tu-Game that reflects the Mathematica convention of 
%              coalition order, that is, with respect to their size.
%  input:
%  v        -- Tu-Game that reflects the Matlab convention of coalition order.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/29/2010        0.1             hme
%   10/27/2012        0.3             hme
%

N=length(v);
[~, n]=log2(N); 
if (2^n-1)~=N, error('Game has not the correct size!'); end
S=1:N;
mg=sortsets(S,n);
w=v(mg);

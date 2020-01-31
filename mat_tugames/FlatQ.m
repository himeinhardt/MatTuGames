function fQ=FlatQ(v,tol)
% FLATQ checks if the game v is flat or not
%
% Usage: eQ=EssentialQ(v) 
%
% Define variables:
%  output:     
%  fQ       -- Returns 1 (ture) whenever the game is flat,
%              otherwise 0 (false)
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  tol      -- A tolerance value, defualt is tol=10^6*eps.
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/20/2014        0.5             hme

if nargin <2
  tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
fQ=abs(v(N))<tol;

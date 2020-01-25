function hsQ=hypergraphQ(hs,n)
% HYPERGRAPHQ returns 1 whenever a coalition structure hs is 
% a hypergraph communication situation.
%
%  Usage: hsQ=hypergraphQ(hs)
%
% Define variables:
%  output:
%  hsQ      -- Returns 1 (true) or 0 (false).
%  input:
%  hs       -- A hypergraph communication situation like
%              [3 6 7 24] for {[1 2],[2 3],[1 2 3], [4 5]}.
%              This returns 1 (true).
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/01/2013        0.4             hme
%   05/15/2014        0.5             hme
%    
    
N=2^n-1;
J=1:n;
pl=2.^(J-1);
hs=sort(hs);
uiQ=any(ismembc(hs,pl));
if uiQ==1
    hsQ=false;
else
    hsQ=all(hs<=N);
end    
    

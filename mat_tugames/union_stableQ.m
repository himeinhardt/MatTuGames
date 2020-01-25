function usQ=union_stableQ(cs)
% UNION_STABLEQ returns 1 whenever a coalition structure cs is 
% union stable.
%
%  Usage: usQ=union_stableQ(cs)
%
% Define variables:
%  output:
%  usQ      -- Returns 1 (true) or 0 (false).
%  input:
%  cs       -- A communication situation like [3 5 6]
%              for {[1,2],[1 3],[2 3]}, this returns 0 (false).
%              Or a union stable structure like [1 7 14 15] for {[1],[1 2 3],
%              [2 3 4],[1 2 3 4]}. This returns 1 (true).
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/26/2013        0.4             hme
%    
    
lcs=length(cs);
k=1;
for ii=1:lcs-1
    for jj=ii+1:lcs
    u=bitand(cs(ii),cs(jj));
     if u > 0
        o=bitor(cs(ii),cs(jj));
        us(k)=any(o==cs);
        k=k+1;
     else
       us(k)=true;
       k=k+1;   
     end
    end
end    
usQ=all(us);

function usQ=p_union_stableQ(cs)
% P_UNION_STABLEQ returns 1 whenever a coalition structure cs is 
% union stable.
%
%  Usage: usQ=p_union_stableQ(cs)
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
bd=lcs-1;
parfor ii=1:bd
    ius=false(1,lcs);
    for jj=1:lcs
    u=bitand(cs(ii),cs(jj));
     if u > 0
        o=bitor(cs(ii),cs(jj));
        ius(jj)=any(o==cs);
     else
        ius(jj)=true;
     end
    end
    us(ii,:)=ius;
end
usQ=all(all(us));

function [bf df]=UnionStableBasis(cs,n)
% UNIONSTABLEBASIS determines a basis of a union stable system
%
% Usage: [bf df]=UnionStableBasis(cs,n)
%
% Define variables:
%  output:
%  bf       -- A basis of a union stable system.
%  df       -- Set that can be written as the union of two 
%              disjoint sets/coalitions.
%
%  input:
%  cs       -- A union stable system like [1 7 14 15] 
%              for {[1],[1 2 3], [2 3 4],[1 2 3 4]} 
%  n        -- Number of the player set.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/13/2015        0.5             hme
%


if nargin < 2
   N=cs(end);
   [~, n]=log2(N);
end

k=1;
J=1:n;
pl=2.^(J-1);
lcs=length(cs);
int=0:-1:1-n;
%if strcmp(str,'us')
 for ii=1:lcs-1  % Constructing the set D(F).  
    for jj=ii+1:lcs
    u=bitand(cs(ii),cs(jj));
     if u > 0
        g=bitor(cs(ii),cs(jj));
        gcs=any(g==cs);
        eQ=(g==[cs(ii),cs(jj)]);
        gQ=eQ==false(1,2);
        if gcs & gQ
          df(k)=g;
          k=k+1;
        end
     end
    end
 end

 % Constructing the basis set B(F).  
df=sort(unique(df));
bf=setdiff(cs,df);
cf=setdiff(bf,pl);

function us=genUnionStable(cs)
% GENUNIONSTABLE generates from a coalition structure cs the 
% union stable system.
%
%  Usage: us=genUnionStable(cs)
%
% Define variables:
%  output:
%  us       -- A union stable system like [3 5 6 7]
%  input:
%  cs       -- A communication situation like [3 5 6]
%              for {[1,2],[1 3],[2 3]}.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/25/2013        0.4             hme
%    
    
lus=length(cs);
if lus>1
us=cs;
  while 1
      cs=SelectingUnions(cs);
      if isempty(cs)
         us=cs;
         break;
      end    
      lc=length(cs);
      lu=length(us);
      if lc==lu
         if us==cs
            break;
         end
      end
      us=cs;
  end
else 
us=cs; 
end        
%--------------------------------
function us=SelectingUnions(cs)    
    
lcs=length(cs);

k=1;    
for ii=1:lcs-1
  for jj=ii+1:lcs
    u=bitand(cs(ii),cs(jj));
    if u > 0
        us{k}=bitor(cs(ii),cs(jj));
        us{k}(2:3)=[cs(ii),cs(jj)];
        k=k+1;
    else
        us{k}=[];
        k=k+1;
    end
  end
end
us=sort(unique(cell2mat(us)));

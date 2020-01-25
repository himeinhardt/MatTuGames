function us=p_genUnionStable(cs)
% P_GENUNIONSTABLE generates from a coalition structure cs the 
% union stable system using Matlab's PCT.
%
%  Usage: us=p_genUnionStable(cs)
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
%   07/26/2013        0.4             hme
%    
    
lus=length(cs);
if lus>1
us=cs;
  while 1
      cs=p_SelectingUnions(cs);
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
function us=p_SelectingUnions(cs)    
    
lcs=length(cs);
us=cell(1,lcs);
bd=lcs-1;
parfor ii=1:bd
  tus=[];
  for jj=1:lcs
    u=bitand(cs(ii),cs(jj));
    if u > 0
        tus=[tus,cs(ii),cs(jj),bitor(cs(ii),cs(jj))];
    end
  end
%  disp(tus)
  us{ii}=tus;
end
us=sort(unique(cell2mat(us)));

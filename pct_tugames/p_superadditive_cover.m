function v_st=p_superadditive_cover(v)
% P_SUPERADDITIVE_COVER computes from game v its superadditive cover using Matlab's PCT.
% It requires Integer Partition Generator from
% 
% SOURCE: http://www.mathworks.com/matlabcentral/fileexchange/36437-integer-partition-generator
%
% It takes some time to finish for n=7. Needs a lot of memory , for instance, 45 GB for n=7.
%
% Usage: v_st=p_superadditive_cover(v)
% Define variables:
%
%  output:
%  v_st     -- The superadditive cover of game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/24/2016        0.8             hme

N=length(v);
[~, n]=log2(N);
v_st=zeros(1,N);
it=0:-1:1-n;
jj=1:n;
parfor S=1:N            % We have to rethink about the partition of integers.
   pt=intpartitions(S); % the partition consumes most of the computing time.
   sm=numel(pt);
   vp=zeros(1,sm);
   idx=[];
   for k=1:sm
       ui=unique(pt{k});
       sQ=sum(ui);
       if sQ==S;
          ps=bitget(S,jj)==1;
          lp=ps*ones(n,1);
          sS=double(ui);
          ci=rem(floor(sS(:)*pow2(it)),2)==1;
          sci=size(ci,1);
          ad=ones(1,sci)*ci;
          if all(ad==ps)
             vp(k)=sum(v(ui));
          else
             idx=[idx,k]; 
          end
       else
          idx=[idx,k]; 
       end
   end
   vp(idx)=[];
   v_st(S)=max(vp);
end

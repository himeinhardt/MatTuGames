function c_v=p_mcst_game(cm)
% P_MCST_GAME computes from a cost matrix the corresponding mcst game 
% using Matlab's PCT. 
%
% Usage: c_v=p_mcst_game(cm)
%
% Define variables:
%  output:
%  c_v      -- A cost game v of length 2^n-1.
%            
%  input:
%
%  cm       -- A square cost matrix (n+1xn+1) derived from a minimum cost
%              spanning tree problem. For instance, for a four
%              person game the size of the matrix must be (5x5). The source
%              is player 1.
% 
%  Example:
%  Let a spanning tree system be represented by the following cost matrix
%  cm = 
%
%     0     2     2     6
%     2     0     1     2
%     2     1     0     2
%     6     2     2     0
%
% then the three person minimum cost spanning tree game is given by
%
%  c_v=p_mcst_game(cm)
%
%  c_v =
%
%     2     2     3     6     4     4     5
%
%  Notice, that the matrix cm must be fully occupied with positive numbers but
%  for the diagonal with zeros.
%
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/27/2016        0.8             hme
%


[m1,m2]=size(cm);
if m1~=m2
   error('Matrix is not square'); 
else
  n=m1;
end    
N=2^n-1;
kk=1:n;
S=1:N;
c_v=zeros(1,N);
for ii=1:n
    mS(:,ii)=bitget(S,ii)==1;
end
parfor sS=1:N
    pl=kk(mS(sS,:));
    lp=length(pl);
    pm=zeros(n);
    for ii=1:lp-1
        for jj=ii+1:lp
            if cm(pl(ii),pl(jj))>0
               pm(pl(ii),pl(jj))=cm(pl(ii),pl(jj));
            end    
        end    
    end
    id=find(pm(1,:)>0);
    if isempty(id)==0
       lid=length(id);

       sw=0;
       for k=1:lid
           id2=pm(1:id(k)-1,id(k))>0;
           cP2=pm(id2,id(k));   
           mx=min(cP2);
           sw=mx+sw;
       end


    else
       sw=inf;
    end

    c_v(sS)=sw;
end
idx=bitget(S,1)==0; % the source must be part of the coalition.
c_v(idx)=[]; % deleting the coalitions which have no connection to
             % the source.
c_v(1)=[];   % deleting the representative of the empty set.

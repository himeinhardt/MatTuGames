function v=weightedGraphGame(E,c)
% WEIGHTEDGRAPHGAME computes from a weighted graph problem (E,c) a
% weighted graph TU game.
%
% Usage: v=weightedGraphGame(E,c)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  E        -- An edge matrix of size (lx2) or a cell of numel l.
%              The source must be given by 1, and the sink by the
%              number of the player set. However, the edge matrix
%              can also be of size (lx3) then c can be empty. 
%  c        -- A weights vector for each edge. Can be set
%              to [] if the edge matrix has size(lx3). The weights
%              vector will then be retrieved from E(:,3).
%
%
% Example:
% Define a matrix of edges given by
% E =
%   1   1   1   1   2   2   3   4
%   2   3   4   5   3   4   4   5
%
% or equivalently
%
% E={[1 2];[1 3];[1 4];[1 5];[2 3];[2,4];[3 4];[4 5]}
% Do not forget to set the separator by semicolon (;) not by comma.
%
% We have 8 edges here. Furthermore, in total we have 5 vertices, and one
% source given by number 1, and a sink by number 5.
% Then specify the weighted vector of each edge.
% w=[2 5 4 3 3 6 2 4];
%
% Here, we have a player set of {1,2,3,4,5}, hence n=5.
%
% Now, invoke
%  v=weightedGraphGame(E,w)
% to get
% v =
%     0    0    2    0    5    3   10    0    4    6   12    2   11   11   22    0    3    0    5    0    8    3   13    4   11   10   19    6   18   15   29
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%  
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/27/2017        0.9             hme
%

if nargin <2
   error('The weights vector must be given');
end
% Reformatting edge matrix E. 
if iscell(E)
   E=cell2mat(E);
% if the source in E is given by zero.
   me=min(E(:,1));
   if me(1)==0
      s=1;
      E=E+s;
      n=max(E(:,2));
   else   
      n=max(E(:,2)); 
   end
else
   n=length(E); 
   [m,cs]=size(E);
   if m==2
      E=E';
      me=min(E(:,1));
      if me(1)==0
         E0=E(1:2,:);
         s=1;
         n=max(E(:,2))+s;
         E=E0+s;
      else   
         n=max(E(:,2));          
      end
   else
      me=min(E(:,1));
      if me(1)==0
         E0=E(:,1:2);
         s=1; 
         n=max(E(:,2))+s;
         E=E0+s;
      else   
         n=max(E(:,2));          
      end
   end
end
% Pre-allocation of variables.
  N=2^n-1;
  k=1:n;
  v=zeros(1,N);
  cl=sum(2.^(E-1),2);
% Getting weighted graph game.  
  for S=1:N
    A=k(bitget(S,k)==1);
    lA=length(A);
    if lA>1
       for kk=1:lA-1
           for jj=kk+1:lA
               sl=2^(A(kk)-1)+2^(A(jj)-1);
               w=c(find(cl==sl));
               if isempty(w)==0
                  v(S)=v(S)+w;
               end   
           end    
       end
    end   
  end    
end    

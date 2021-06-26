function [v,mW,wC,n,E]=NetworkMajorityGame(E,th,c)
% NETWORKMAJORITYGAME computes from a network problem (E,c,th) a
% network majority TU game (simple game).
%
%
% Source: M. Holler and F. Rupp, Power in Networks: The Medici,2020.
%
%
% Usage: [v,mW]=NetworkMajorityGame(E,c,th)
%
% Define variables:
% output:
%  v        -- A weighted majority Tu-Game v of length 2^n-1. 
%  mW       -- minimal winning coalitions.
%  wC       -- set of winning coalitons.
%  n        -- number of players involved.
%  E        -- Collection of edges represented as an edge matrix.
%
% input: 
%  E        -- An edge matrix of size (lx2) or a cell of numel l.
%              The source must be given by 1, and the sink by the
%              number of the player set. However, the edge matrix
%              can also be of size (lx3) then c can be empty.
%  th       -- Threshold to pass a bill, 1<th<n. 
%  c        -- A weighted vector for each player. 
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
% Set th=3 for the threshold to pass a bill.
% Optinally, one can specify the voting weight of each player by a weighted vector, for instance,
% by w=[2 5 6 2 4]. Otherwise, the default of equal weight is assumed.
%
% Here, we have a player set of {1,2,3,4,5}, hence n=5.
%
% Now, invoke
%  v=NetworkMajorityGame(E,th)
% to get
% v =
%     0   0   0   0   0   0   1   0   0   0   1   0   1   1   1   0   0   0   1   0   1   0   1   0   1   1   1   1   1   1   1
% 
% For the example with weight vector w and th=10, we get
%  [v,mW]=NetworkMajorityGame(E,10,w)
%
% v =
%     0   0   0   0   0   1   1   0   0   0   0   0   1   1   1   0   0   0   1   0   1   1   1   0   0   1   1   1   1   1   1
%
% mW =
%   6    13    19    21    26    28
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/30/2020        1.9             hme
%

if nargin <2
   error('The threshold (quorum) must be given');
end
if iscell(E)
   E=cell2mat(E);
   rw=size(E);
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
   [m,rw]=size(E);
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
N=2^n-1;
k=1:n;
if nargin<3
   %w=ones(1,n);
   uq=unique(E);
   w=ismember(k,uq);   
 else
   w=c;
end
% Computing the measure.
Xm=w(1); for ii=2:n, Xm=[Xm w(ii) Xm+w(ii)]; end
v=false(1,N);
%% Assigning the coalitonal values.
%% No self-loops are allowed.
for S=1:N
    A=k(bitget(S,k)==1);
    lA=length(A);
    if lA>1 && Xm(S)>=th
       it=1;
       q=[];
       for jj=1:lA-1
	  for ii=jj+1:lA
              for kk=1:rw
                  if all([A(jj),A(ii)]==E(kk,:))
                     q(:,it)=[A(jj);A(ii)];
		     it=it+1;
                     break
	          end
              end
          end
       end
       %% Checking if the selected subgraph 
       %% is feasible and connected.
       if isempty(q)==0
          g=graph(q(1,:),q(2,:));
          bins = conncomp(g);
          lb=length(bins);
          pq=unique(q(:)');
          lq=length(pq);
          pl=1:lb;
          if lq==lA
             sbg=pl(bins==A(1));
             if all(ismember(A,sbg))==1
	        v(S)=true;
             end
          end
       end	  
    end
end    
S=1:N;
% pre-selected winning coalitions.
wC=S(v);
lwC=length(wC);
mWS=zeros(1,lwC);
%% Minimal winning coalitions.
for k=1:lwC
  sS=SubSets(wC(k),n);
  [vl, idx]=max(v(sS));
  smS=sS(idx);
  if smS==wC(k)
      mWS(k)=wC(k);
  end
end

mW=mWS(logical(mWS>0));
% winning coalitions.
wC=winning_coalitions(mW,n);
v=ismembc(S,wC);

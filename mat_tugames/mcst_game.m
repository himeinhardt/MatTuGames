function CV=mcst_game(cm,method)
% MCST_GAME computes from a cost matrix the corresponding mcst game. 
% Using Prim's or Kruskal’s algorithm.
%
% Usage: CV=mcst_game(cm)
%
% Define field variables:
%  output:
%  c_v      -- A cost game v of length 2^n-1.
%  tslm     -- Matrix of minimum cost tree solutions.
%  trm      -- Cell of tree solution for each m.c.s.t. problem (C_S,S_0).
%            
%  input:
%
%  cm       -- A square cost matrix (n+1xn+1) derived from a cost
%              spanning graph problem. For instance, for a four
%              person game the size of the matrix must be (5x5). The source
%              is player 1.
%  method   -- A string to define the method. Admissible methods are:
%              'prim' 
%              'kruskal'
%              The default method is 'prim'.
%
% 
%  Example:
%  Let a spanning graph system be represented by the following cost matrix
%  cm = 
%
%     0     2     2     6
%     2     0     1     2
%     2     1     0     2
%     6     2     2     0
%
% then the three person minimum cost spanning tree game is given by
%
%  CV=mcst_game(cm)
%
%  CV.c_v =
%
%     2     2     3     6     4     4     5
%
%  Notice, that the matrix cm must be fully occupied with non-zeros but
%  for the diagonal. It represents a complete undirected graph.
%
%
  

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/27/2021        1.9             hme
%

if nargin < 2
   method='prim';
end	

[m1,m2]=size(cm);
if m1~=m2
   error('Matrix is not square'); 
else
  n=m1;
end    
N=2^n-1;
n1=n-1;
N1=2^n1-1;
kk=1:n;
S=1:N;
c_v=zeros(1,N);
S1=1;
tslm=zeros(N1,n1);
for ii=1:n
    mS(:,ii)=bitget(S,ii)==1;
end
upe=true(n);
ann=n^2;
ar=1:ann;
szA=[n,n];
A=reshape(ar,szA);
UA=triu(upe,1);
ind=A(UA)';
trm=cell(N1,1);

for sS=1:N
    pl=kk(mS(sS,:));
    lp=length(pl);
    pm=zeros(n);
    %% Setting up cost matrix for coalition sS.
    for ii=1:lp-1
        for jj=ii+1:lp
            if cm(pl(ii),pl(jj))>0
               pm(pl(ii),pl(jj))=cm(pl(ii),pl(jj));
            end    
        end    
    end
    %% Finding players that connect to the source 
    %% to cheapest costs.
    zr=pm(1,:);
    id=find(zr>0);
    if isempty(id)
       continue;
    end
    czr=zr(zr>0);
    mc=min(czr);
    mid=id(czr==mc);
    ccl=czr(czr==mc);
    zidx=[];
    lmd=length(mid);
    %% Checking if there is a cheaper indirect connection. 
    if lmd>1
      for ii=1:lmd
          zcl=pm(:,mid(ii));
          czcl=zcl(zcl>0);
          mzcl=min(czcl);
          if mzcl<mc
             zidx=[zidx,ii]; 
          end    
      end
      if isempty(zidx)==0
         mid(zidx)=[];
         ccl(zidx)=[];
      end    
    end
   %% Finding minimum cost tree solution.
    if isempty(id)==0    
       wghs=pm(UA)';
       slc=wghs>0;
       wghs=wghs(slc);
       [s,t]=ind2sub(szA,ind);
       s=s(slc);
       t=t(slc);
       G = graph(s,t,wghs);
       if strcmp(method,'prim')
       %% Prim’s algorithm.
          [T,pred] = minspantree(G,'Method','dense');
       elseif strcmp(method,'kruskal')	  
       %% Kruskal’s algorithm.
          [T,pred] = minspantree(G,'Method','sparse');
       else
          [T,pred] = minspantree(G,'Method','dense');
       end	       
       wghs2=T.Edges.Weight';
       sw=sum(wghs2);
       %% Formatting output
       mE=T.Edges.EndNodes;
       mE2=mE'-1;
       trm{sS}=mE;
       la=size(mE);
       t1=pred-1;
       t1(1)=[];
       t1=t1(~isnan(t1));
       mt=[1:la(1);t1];
       for ll=1:la(1)
	  tk=mt(:,ll);
	  stk=sort(tk);
          for jj=1:la(1)
              if all(tk==mE2(:,jj))
		 tslm(S1,ll)=wghs2(jj);     
		 continue;     
              elseif all(stk==mE2(:,jj))
	         tslm(S1,ll)=wghs2(jj); 	      
	         continue;
              end		 
	  end	  
       end
       S1=S1+1;   
    else
       sw=inf;
    end
    c_v(sS)=sw;
end
idx=bitget(S,1)==0; % the source must be part of the coalition.
c_v(idx)=[]; % deleting the coalitions which have no connection to
%             % the source.
c_v(1)=[];   % deleting the representative of the empty set.
CV.c_v=c_v;
CV.tslm=tslm;
CV.trm=trm(~cellfun('isempty',trm));

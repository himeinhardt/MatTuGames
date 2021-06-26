function ICM=IrredCostMatrix(cm,mE)
% IRREDCOSTMATRIX computes from a cost matrix and a solution tree the irreducible cost matrix. 
%
% Usage: icm=IrredCostMatrix(cm,mE)
% 
%
%  Source: D. Granot and G. Huberman, Minimum Cost Spanning Tree Games, Mathematical Programming 21 (1981), 1-18.
%    
% Define structure field variables:
%  output:
%  icm      -- An irreducible square cost matrix (n+1xn+1) derived from cost matrix cm and
%              and a minimal tree solution tr.     
%  cm       -- The default cost matrix (input).
%            
%  input:
%
%  cm       -- A square cost matrix (n+1xn+1) derived from a minimum cost
%              spanning tree problem. For instance, for a four
%              person game the size of the matrix must be (5x5). The source
%              is player 1.
%  mE       -- Tree solution of a minimum cost spanning tree problem 
%              in matrix form (edge matrix).
%
%    
%  Example:
%  Let a spanning tree system be represented by the following cost matrix
% cm =
%
%     0     2     4     7     2     3
%     2     0     3     1     3     4
%     4     3     0     4     2     3
%     7     1     4     0     6     5
%     2     3     2     6     0     2
%     3     4     3     5     2     0
%
% and a minimal tree solution of the m.c.s.t. problem is given by
%
% mE = 
%
%     1     1     2     3     5
%     2     5     4     5     6
%
% then the irreducible cost matrix is obtained by
%
% ICM=IrredCostMatrix(cm,mE);
%
% ICM.icm 
%
% ans =
%
%     0     2     2     2     2     2
%     2     0     2     1     2     2
%     2     2     0     2     2     2
%     2     1     2     0     2     2
%     2     2     2     2     0     2
%     2     2     2     2     2     0
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
    

[sz1,sz2]=size(cm);
icm=zeros(sz1,sz1);
if sz1==sz2
   n=sz1;
else
   error('Cost Matrix is not symmetric');	
   return	
end	

N=2^n-1;
n1=n-1;
N1=2^n1-1;

la=size(mE);
if la(2)==2
   mE=mE';
end	
zr=cm(1,:);
n1=max(mE(:));
la=size(mE);

for kk=1:la(2)
   idx=mE(1,kk);
   idy=mE(2,kk); 
   if icm(idx,idy)==0
       cv(kk)=cm(idx,idy);
       icm(idx,idy)=cm(idx,idy);
       icm(idy,idx)=cm(idy,idx);
   end
end

upe=true(n);
ann=n^2;
ar=1:ann;
szA=[n,n];
A=reshape(ar,szA);
UA=triu(upe,1);
ind=A(UA)';
wghs=icm(UA)';
slc=wghs>0;
wghs=wghs(slc);
[s,t]=ind2sub(szA,ind);
s=s(slc);
t=t(slc);
G = graph(s,t,wghs);

for ii=1:n-1;
    for jj=ii+1:n
        if icm(ii,jj)==0	    
           [P,d,edgepath] = shortestpath(G,ii,jj);
           shw=G.Edges(edgepath,:);
           icm(ii,jj)=max(shw.Weight);
           icm(jj,ii)=max(shw.Weight);
        end
    end	
end    
ICM.icm=icm;
ICM.cm=cm;

    
    

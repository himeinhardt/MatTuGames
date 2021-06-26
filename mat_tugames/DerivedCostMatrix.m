function DCM=DerivedCostMatrix(cm,cs)
% DERIVEDCOSTMATRIX computes from a cost matrix and a partition of the player set N the corresponding derived cost matrix. 
%
% Usage: dcm=DerivedCostMatrix(cm,cs)
% 
%
%  Source: D. Granot and G. Huberman, Minimum Cost Spanning Tree Games, Mathematical Programming 21 (1981), 1-18.
%    
% Define structure field variables:
%  output:
%  dcm      -- A derived square cost matrix (n+1xn+1) derived from cost matrix cm and
%              and partition cs.     
%  icm      -- Cell array of induced cost matrices from (dcm,cs). 
%            
%  input:
%
%  cm       -- A square cost matrix (n+1xn+1) derived from a minimum cost
%              spanning tree problem. For instance, for a four
%              person game the size of the matrix must be (5x5). The source
%              is player 1.
%  cs       -- A coalition structure provided as partition of N like [1 6].
%              Restriction: The length of cs cannot be larger than 2.    
%
%    
%  Example:
%  Let a spanning tree system be represented by the following cost matrix
%  cm = 
%
%   0   2   4   7   2   3
%   2   0   3   1   3   4
%   4   3   0   4   2   3
%   7   1   4   0   6   5
%   2   3   2   6   0   2
%   3   4   3   5   2   0
%
% and the partition cs of the player set is given by
%
% cs=[5 26]
%
% then the derived cost matrix is obtained by
%
% dcm=DerivedCostMatrix(cm,cs)
%
%   0   2   3   4   2   3
%   2   0   3   1   3   4
%   3   3   0   4   2   3
%   4   1   4   0   6   5
%   2   3   2   6   0   2
%   3   4   3   5   2   0    
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
    
N=bitor(cs(1),cs(2));
[~, n]=log2(N);
it=0:-1:1-n;    
slcP=rem(floor(cs(:)*pow2(it)),2)==1;
J=1:n;
[rw,cl]=size(slcP);
sz=slcP*ones(n,1);
dcm=cm;
icm=cell(1,rw);

for kk=1:rw
   sP=J(slcP(kk,:));
   n2=setdiff(J,sP);
   idx=sP+1;
   idy=n2+1;
   for jj=1:sz(kk)
       c0=cm(idx(jj),1);
       mc1=min(cm(idx(jj),idy));
       dcm(idx(jj),1)=min([c0,mc1]);
       dcm(1,idx(jj))=min([c0,mc1]);
   end
   eidx=[1,idx];
   icm{kk}=dcm(eidx,eidx);
end
DCM.dcm=dcm;
DCM.icm=icm;


    
    

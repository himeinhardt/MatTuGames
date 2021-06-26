function MTR=MTRCostMatrix(cm,mE)
% MTRCOSTMATRIX computes from a cost matrix and a solution tree the cost matrix of a minimal spanning tree. 
%
% Usage: MTR=MTRCostMatrix(cm,mE)
% 
%
%  Source: D. Granot and G. Huberman, Minimum Cost Spanning Tree Games, Mathematical Programming 21 (1981), 1-18.
%    
% Define structure field variables:
%  output:
%  mtr      -- A square cost matrix (n+1xn+1) derived from cost matrix cm and
%              and a minimal tree solution mE.     
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
% cm
%
% cm =
%      0     2     4     7     2     3
%      2     0     3     1     3     4
%      4     3     0     4     2     3
%      7     1     4     0     6     5
%      2     3     2     6     0     2
%      3     4     3     5     2     0
%
% and a tree solution of m.c.s.t. problem given by
%
% mE =
%      1     1     2     3     5
%      2     5     4     5     6
%
% then the  cost matrix is obtained by
%
% MTR=MTRCostMatrix(cm,Em);
% MTR.mtr
%
% ans =
%     0     2     0     0     2     0
%     2     0     0     1     0     0
%     0     0     0     0     2     0
%     0     1     0     0     0     0
%     2     0     2     0     0     2
%     0     0     0     0     2     0
%
%  This matrix can then be used to plot the minimal tree solution of a m.c.s.t. problem
%  with the function PlotCostGraph(MTR.mtr).
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
MTR.mtr=icm;
MTR.cm=cm;

    
    

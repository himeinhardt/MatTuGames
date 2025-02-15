function PKM=cs_GetPrk(v,tol)
% CS_GETPRK computes a pre-kernel element from each possible
% partition of N.
%
% Usage: PKM=cs_GetPrk(v)
%
% Define variables:
%  Structure element:
%  pkm      -- Matrix of Pre-Kernel elements (output). Ordered by 
%              the cell of partitions clm.
%  pkQ      -- Returns a list of ones (true) and zeros (false).
%              Each component indicates whether the solution is
%              a pre-kernel element of underlying coalition structure.
%  clm      -- List of Partitions of N.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/01/2017        0.9             hme
%                


if nargin<2
   tol=10^7*eps;
end

N=length(v);
[~, n]=log2(N);
clm=GetPartitions(N);
sc=numel(clm); %% Is a number from the Bell triangle.
pkm=zeros(sc,n);
pkQ=false(1,sc);
for k=1:sc
    pkm(k,:)=cs_PreKernel(v,clm{k});
    pkQ(k)=cs_PrekernelQ(v,clm{k},pkm(k,:),tol);
end    

PKM.pkm=pkm;
PKM.pkQ=pkQ;
PKM.clm=clm;

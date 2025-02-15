function PKM=cs_GetPrk(clv,tol)
% CS_GETPRK computes a pre-kernel element from each possible
% partition of N.
%
% Usage: PKM=clv.cs_GetPrk()
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
%  clv      -- TuGame class object.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/22/2017        0.9             hme
%                


if nargin<2
   tol=10^7*eps;
end

N=clv.tusize;
n=clv.tuplayers;
clm=GetPartitions(N);
sc=numel(clm); %% Is a number from the Bell triangle.
pkm=zeros(sc,n);
pkQ=false(1,sc);
for k=1:sc
    pkm(k,:)=clv.cs_PreKernel(clm{k});
    pkQ(k)=clv.cs_PrekernelQ(clm{k},pkm(k,:),tol);
end    

PKM.pkm=pkm;
PKM.pkQ=pkQ;
PKM.clm=clm;

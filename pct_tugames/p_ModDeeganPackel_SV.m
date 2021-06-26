function dpidx=p_ModDeeganPackel_SV(sv)
% P_ModDeeganPackel_SV computes the Deegan-Packel index from a simple game to construct the set of minimal winning coalitions
% using Matlab's PCT.
%
% Usage: dpidx=p_ModDeeganPackel_SV(sv)
% Define variables:
%  output:
%  dpidx    -- The modified Deegan-Packel index.
%
%  input:
%  sv       -- A simple game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/01/2020        1.9             hme
%

N=length(sv);
[~, n]=log2(N);
sWCk=zeros(1,n);
[mWC,WC]=p_getMinimalWinning(sv);
m=length(WC);
dpidx=zeros(1,n);
it=0:-1:1-n;
mat=rem(floor(mWC(:)*pow2(it)),2);
parfor k=1:n;
    if any(mat(:,k)==1)
       mWCk=WC(bitget(WC,k)==1);
       sk=length(mWCk);
       A=zeros(sk,1);
       for jj=1:n, 
           A(:,jj) = bitget(mWCk,jj);
       end
       ssk=sum(1./(A*ones(n,1)));
       dpidx(k)=ssk/m;
   end    
end


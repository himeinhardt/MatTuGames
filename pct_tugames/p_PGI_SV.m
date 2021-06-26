function pgi=p_PGI_SV(sv)
% P_PGI_SV computes the public good index from a simple game to determine the set of minimal winning coalitions
% using Matlab's PCT.
%
% Usage: gpi=p_PGI_SV(sv)
% Define variables:
%  output:
%  pgi      -- The Holler/Public good index.
%
%  input:
%  sv       -- Data array of a single game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/03/2020        1.9             hme
%


N=length(sv);
[~, n]=log2(N);
sWCk=zeros(1,n);
mWC=p_getMinimalWinning(sv);
parfor k=1:n;
    mWCk=mWC(bitget(mWC,k)==1);
    sWCk(k)=length(mWCk);
end
pgi=sWCk./(sWCk*ones(n,1));


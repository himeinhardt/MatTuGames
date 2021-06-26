function mpgi=p_ModPGI_SV(clv)
% P_ModPGI_SV computes the modified public good index from a simple game to determine the set of minimal winning coalitions
% using Matlab's PCT.
%
% Usage: mgpi=clv.p_ModPGI_SV()
% Define variables:
%  output:
%  mpgi     -- The modified Holler/Public good index.
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
%   10/03/2020        1.9             hme
%


N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;
if strcmp(gt,'sv')
else
   error('Wrong game type!. Game must be a simple game!')
end
sWCk=zeros(1,n);
[mWC,wC]=clv.p_getMinimalWinning();
it=0:-1:1-n;
mat=rem(floor(mWC(:)*pow2(it)),2);
parfor k=1:n;
    if any(mat(:,k)==1)
       mWCk=wC(bitget(wC,k)==1);
       sWCk(k)=length(mWCk);
    end   
end
mpgi=sWCk./(sWCk*ones(n,1));


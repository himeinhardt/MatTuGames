function pgi=PGI_SV(clv)
% PGI_SV computes the public good index from a simple game to determine the set of minimal winning coalitions.
%
% Usage: gpi=clv.PGI_SV()
% Define variables:
%  output:
%  pgi      -- The Holler/Public good index.
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
mW=clv.getMinimalWinning();
for k=1:n;
    mWCk=mW(bitget(mW,k)==1);
    sWCk(k)=length(mWCk);
end
pgi=sWCk./(sWCk*ones(n,1));


function bcQ=p_B0_balancedCollectionQ(v,x,tol)
% P_B0_BALANCEDCOLLECTIONQ verifies whether the set of induced
% coalitions is a B0_balanced collection. Checking Kohlberg's criterion.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
% Uses now Dual-Simplex (Matlab R2015a).
%
%
% Usage: bcQ=p_B0_balancedCollectionQ(v,x,tol)
% Define variables:
%  output:
%  bcQ      -- Returns 1 (true) or 0 (false).
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  x        -- A (pre-)imputation of length n.
%  tol      -- Tolerance value. Its default value is set to 10^4*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/06/2017        0.9             hme
%


if nargin < 2
   tol=10^4*eps;
elseif nargin < 3
   tol=10^4*eps;
end

bcQ=false;
N=length(v);
[~, n]=log2(N);
effQ=abs(sum(x)-v(N))<tol;

if effQ==0
   return;
end

exc=excess(v,x);
k=1:n;
ic=2.^(k-1);
if any(exc(ic)>0)
   return;
else
  iex=exc(ic)==0;
end
b0=ic(iex);
pws=PowerSet(b0);
lb0=length(pws);
leb=lb0+1;
cS=cell(1,leb);

exc(N)=[];
[sx,idx]=sort(exc,'descend');
mx=max(sx);
slc=sx>=mx-tol;
S=idx(slc);
parfor k=1:lb0
   cS{k}=unique([S,pws{k}]);
end
cS{leb}=S;



p = gcp(); % get the current parallel pool


for kk=1:leb;
    F(kk) = parfeval(p,@p_balancedSetQ,1,exc,n,cS{kk},tol);
end

%% Borrowed from Matlab Documentation
for kk = 1:leb
    [widx, thisResult] = fetchNext(F);
    if thisResult == 1
        %widx
        bcQ = thisResult;
        % Have all the results needed, so break
        break;
    end
end
% With required result, cancel any remaining futures
cancel(F);



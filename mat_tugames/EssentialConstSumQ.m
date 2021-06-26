function creQ=EssentialConstSumQ(v,tol)
% ESSENTIALCONSTSUMQ checks if v is an essential constanst-sum game, then
% the core is empty.
%
% Usage: creQ=EssenticalConstSumQ(v) 
%
% Define variables:
%  output:     
%  creQ     -- Returns 1 (ture) whenever v is an essential constant-sum game,
%              otherwise 0 (false)
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  tol      -- A tolerance value, defualt is tol=10^6*eps.
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/21/2014        0.5             hme

if nargin <2
  tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N-1;
CN=N-S;
cv=v(CN);
cv(N)=0;
vN=ones(1,N)*v(N);
csQ=all(abs(v+cv-vN)<tol);
eQ=EssentialQ(v);
creQ=eQ & csQ;

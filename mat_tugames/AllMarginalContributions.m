function Mgc=AllMarginalContributions(v)
% ALLMARGINALCONTRIBUTIONS computes all marginal worth vectors 
% of a TU-game v. For n=12, the computation requires at least 360 GB.
%
% Usage: Mgc=AllMarginalContributions(v)
% Define variables:
%  output:
%  Mgc      -- The matrix of marginal contributions.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/29/2010        0.1 beta        hme
%   10/27/2012        0.3             hme
%   09/07/2020        1.9             hme
%   06/08/2023        1.9.1           hme
%                
narginchk(1,1); % check for legal number of input arguments.
%
% Section for assigning variables and pre-allocations.
N=length(v);
[~, n]=log2(N);
pl=1:n;
ci=2.^(pl-1);
spm=perms(ci);
N1=N+1; %% defining empty set.
v1=[v,0];
sz=size(spm);
A=triu(ones(n));
Mgc=zeros(sz);
P=spm*A;
[~,ix]=sort(spm,2);
clear spm;
%% Reordering matrix P w.r.t. ix.
clm=repmat(1:sz(1),1,n);
rw=ix(:)';
idx=sub2ind([n sz(1)],rw,clm);
Pt=P';
P=reshape(Pt(idx),sz(1),n);
clear ix idx rw clm Pt;
% Section marginal worth vectors.
rP1=zeros(sz);
% Initialization of the baseline matrix of permutations.  
for k=1:n
    rP1(:,k)=P(:,k)-ci(k);
end
%% Assigning the empty set.
rP1(rP1==0)=N1;
%% Computing Marginal Contributions.
Mgc=v(P)-v1(rP1);

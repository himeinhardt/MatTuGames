function [Mgc mgc P ix]=p_AllMarginalContributions(v)
% P_ALLMARGINALCONTRIBUTIONS computes all marginal worth vectors 
% of a TU-game v. Using Matlab's PCT. For n=12, one needs
% at least 400 GB.
%
% Usage: [Mgc mgc P ix]=p_AllMarginalContributions(v)
% Define variables:
%  output:
%  Mgc      -- The matrix of marginal contributions.
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
%   05/22/2011        0.1 alpha        hme
%   10/27/2012        0.3              hme
%   08/31/2020        1.9              hme
%                
msg=nargchk(1,1,nargin);
error(msg);

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
rP1=P;
%% requires too much momery for n=12; 
%% negative impact on the performance.
%spmd
% codistributed(P);
% codistributed(rP1);
%end

parfor k=1:n
    rP1(:,k)=P(:,k)-ci(k);
end
%rP1=gather(rP1);
%% Assigning the empty set.
rP1(rP1==0)=N1;
%% Computing Marginal Contributions.
Mgc=v(P)-v1(rP1);

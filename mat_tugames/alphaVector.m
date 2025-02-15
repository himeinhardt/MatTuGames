function [alpv, v]=alphaVector(v,pt,tol)
% ALPHAvECTOR computes recursively an alpha vector from the positive core of game v 
% using Matlab's Optimization Toolbox. Uses dual-simplex (Matlab R2015a).
%
%
% Usage: [alpv, v]=alphaVector(v,pt,tol)
% Define variables:
%  output:
%  alpv     -- Alpha vector computed from the positive core.
%  v        -- A derived Tu-Game v of length 2^n-1.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  pt       -- A partition of N by their unique integers representation.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/02/2015        0.6             hme
%   03/30/2015        0.7             hme
%   04/22/2024        1.9.2           hme
%                


if nargin<3
 tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
% Determining the lower bounds. 
S=1:N;
%k=1:n;
%vi=v(bitset(0,k));
%lb=vi';
lb=[];
% Determining the pre-nucleolus
try
   x=cplex_prenucl(v);
catch
   x=PreNucl(v);
end
exv=max(excess(v,x),0);
ad=additive_game(x);
cex=S(exv>tol);
v(cex)=ad(cex);
[crv,crst]=CddPositiveCoreVertices(v);
Ppcr=Polyhedron(crv);
Ppcr.computeHRep;


% Constructing the LP
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(N,:);
A2=sparse(A1);
lpt=length(pt);
alpv=zeros(1,lpt);
cmat=-A1(pt,:);
ub=[];
%opts=optimset('Simplex','on','LargeScale','off','MaxIter',256);
opts.Simplex='on';
opts.Display='off';
%opts.ActiveSet='on'; % wrong results!
opts.LargeScale='on';
mth1=verLessThan('matlab','24.1.0');
if mth1==0,
    opts.Algorithm='dual-simplex-highs';
else
    opts.Algorithm='dual-simplex';
end
opts.TolFun=1e-10;
opts.TolX=1e-10;
opts.TolRLPFun=1e-10;
% dual-simplex options
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);


for ii=1:lpt
    B1=[-v';v(N)];
    C=cmat(ii,:);
    try
        [x,fmin,exitflag,output,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
    catch
        [x,fmin,exitflag,output,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
    end	
    if exitflag<0
       exitflag
       error('Optimization terminated not successfully');
    end
    pcQ=Ppcr.contains(x);
    if pcQ==0
       warning('Warn:Pcr','Solution found does not belong to the positive Core!');
    end
    alpv(ii)=cmat(ii,:)*x;
    v(pt(ii))=alpv(ii);
end

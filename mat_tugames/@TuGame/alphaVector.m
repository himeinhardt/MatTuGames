function [alpv, v]=cplex_alphaVector(clv,pt,tol)
% ALPHAvECTOR computes recursively an alpha vector from the positive core of game v 
% using Matlab's Optimization Toolbox. Uses dual-simplex (Matlab R2015a).
%
%
% Usage: [alpv, v]=clv.cplex_alphaVector(pt,tol)
% Define variables:
%  output:
%  alpv     -- Alpha vector computed from the positive core.
%  v        -- A derived Tu-Game v of length 2^n-1.
%
%  input:
%  clv      -- TuGame class object.
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
%   02/24/2018        0.9             hme
%   05/25/2024        1.9.2           hme
%                


if nargin<3
 tol=10^6*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
% Determining the lower bounds. 
S=1:N;
%k=1:n;
%vi=v(bitset(0,k));
%lb=vi';
lb=[];
% Determining the pre-nucleolus
try
   x=clv.cplex_prenucl();
catch
   x=clv.PreNucl();
end
exv=max(clv.excess(x),0);
ad=additive_game(x);
cex=S(exv>tol);
v(cex)=ad(cex);
[crv,crst]=clv.CddPositiveCoreVertices();
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
x0=[];
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


%warning('off','all');
%mtv=verLessThan('matlab','9.1.0');
%if mtv==1
%  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
%else
%  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
%end
%warning('on','all');
%options

for ii=1:lpt
    B1=[-v';v(N)];
    C=cmat(ii,:);
    try
        [x,fmin,exitflag,output,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
    catch
        [x,fmin,exitflag,output,lambda]=linprog(C,A2,B1,[],[],lb,ub,x0,opts);
    end    
%    [x,fmin,exitflag,output,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
    if exitflag<0
       exitflag
       error('Optimization terminated not successfully')
    end
    pcQ=Ppcr.contains(x);
    if pcQ==0
       warning('Warn:Pcr','Solution found does not belong to the positive Core!');
    end    
    alpv(ii)=cmat(ii,:)*x;
    v(pt(ii))=alpv(ii);
end


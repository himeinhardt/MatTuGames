function [alpv, v]=gurobi_alphaVector(clv,pt,tol)
% GUROBI_ALPHAvECTOR computes recursively an alpha vector from the positive core of game v 
% using gurobimex.
%
% GUROPI-SOLVER: http://www.gurobi.com
%
% Usage: [alpv, v]=clv.gurobi_alphaVector(pt,tol)
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
%                


if nargin<3
 tol=10^6*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
% Determining the pre-nucleolus
S=1:N;
try
   x=clv.gurobi_prenucl();
catch
   x=clv.PreNucl();
end
exv=max(clv.excess(x),0);
ad=additive_game(x);
cex=S(exv>tol);
v(cex)=ad(cex);
%k=1:n;
%vi=v(bitset(0,k));
%lb=vi';
lb=-inf(n,1);

% Constructing the LP
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(N,:);
A2=sparse(A1);

lpt=length(pt);
alpv=zeros(1,lpt);
cmat=-A1(pt,:);
ub=inf(n,1);

s0='<';
ctype=repmat(s0,1,N);
ctype(N+1)='=';
params.outputflag = 0;
% params.method= 0; % Use primal simplex method.
% params.method= 1; % Use dual simplex method.
% params.method= 2; % Use barrier method.
params.method= 3; % Use concurrent.
% params.method= 4; % Use deterministic concurrent.
params.FeasibilityTol = 1e-9;
params.OptimalityTol = 1e-9;
params.BarConvTol = 1e-10;

model.sense = ctype;
model.vtype = 'C';
model.modelsense = 'min';
model.lb=lb;
model.ub=ub;



for ii=1:lpt
    B1=[-v';v(N)];
    C=cmat(ii,:);
    model.rhs = B1;
    model.A=A2;
    model.obj=C;
    result = gurobi(model,params);
    if strcmp(result.status, 'OPTIMAL')
      exitflag = 1;
    elseif strcmp(result.status, 'ITERATION_LIMIT')
      exitflag = 0;
    elseif strcmp(result.status, 'INF_OR_UNBD')
       params.dualreductions = 0;
       result = gurobi(model, params);
       if strcmp(result.status, 'INFEASIBLE')
          exitflag = -2;
       elseif strcmp(result.status, 'UNBOUNDED')
          exitflag = -3;
       else
          exitflag = nan;
       end
    elseif strcmp(result.status, 'INFEASIBLE')
       exitflag = -2;
    elseif strcmp(result.status, 'UNBOUNDED')
       exitflag = -3;
    else
       exitflag = nan;
    end
    if exitflag<0
       exitflag
       error('Optimization terminated not successfully');
    end
%    x=result.x
    alpv(ii)=result.objval;
    v(pt(ii))=alpv(ii);
    ctype(pt(ii))='=';
end

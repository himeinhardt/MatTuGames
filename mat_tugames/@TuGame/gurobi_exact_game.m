function [v_e xm]=gurobi_exact_game(clv,tol);
% GUROBI_EXACT_GAME computes the exact game from game v using gurobimex.
%
% GUROPI-SOLVER: http://www.gurobi.com 
% 
%  Usage: [v_e xm status]=clv.gurobi_exact_game(tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
%    
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to -10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/25/2016        0.8             hme
%                
if nargin<2
 tol=10^6*eps;
end

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

if clv.CddCoreQ(tol)==0
  error('Core is empty!');
 else
end

k=1:n;
vi=v(bitset(0,k));

S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat;
A2=sparse(A1);
B1=-v';
ub=inf(1,n);
lb=vi';
v_e=zeros(1,N);
xm=zeros(N,n);
s0='<';
ctype=repmat(s0,1,N);
ctype(N)='=';
params.outputflag = 0;
% params.method= 0; % Use primal simplex method.
% params.method= 1; % Use dual simplex method.
params.method= 2; % Use barrier method.
% params.method= 3; % Use concurrent.
% params.method= 4; % Use deterministic concurrent.
params.TimeLimit = 1000;

A2=sparse(A1);
model.A=A2;
model.rhs = B1;
model.sense = ctype;
model.vtype = 'C';
model.modelsense = 'min';
model.lb=lb';
model.ub=ub';


for S=1:N
    C=PlyMat(S,:);
    model.obj=C;
    result = gurobi(model,params);
    x=result.x;
    v_e(S)=C*x;
    xm(S,:)=x;
end    



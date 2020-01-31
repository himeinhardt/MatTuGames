function [v_e xm]=p_gurobi_exact_game(v,tol);
% P_GUROBI_EXACT_GAME computes the exact game from game v using gurobimex
% and Matlab's PCT.
%
% GUROPI-SOLVER: http://www.gurobi.com 
% 
%  Usage: [v_e xm status]=p_gurobi_exact_game(v,tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
%    
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to -10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/25/2016        0.8             hme
%   05/16/2019        1.1             hme
%                
if nargin<2
 tol=10^6*eps;
end

N=length(v);
[~,n]=log2(N);

if CddCoreQ(v,tol)==0
  error('Core is empty!');
 else
end

k=1:n;
vi=v(bitset(0,k));

S=1:N;
parfor k=1:n, PlyMat(:,k) = bitget(S,k);end
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
params.Threads = 1;
C=zeros(n,1); % we need to specify the whole model.
parfor kk=1:N
  A2=sparse(A1);
  model(kk).A=A2;
  model(kk).rhs = B1;
  model(kk).sense = ctype;
  model(kk).vtype = 'C';
  model(kk).modelsense = 'min';
  model(kk).lb=lb';
  model(kk).ub=ub';
  model(kk).obj=C;
end
clear S;
parfor S=1:N
    C=PlyMat(S,:);
    model(S).obj=C';
    result = gurobi(model(S),params);
    xx=result.x;
    v_e(S)=C*xx;
    xm(S,:)=xx;
end    



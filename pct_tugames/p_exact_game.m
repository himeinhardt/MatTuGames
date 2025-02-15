function [v_e xm status]=p_exact_game(v,tol);
% P_EXACT_GAME computes the exact game from game v using 
% the Matlab's Optimization toolbox and PCT. Uses Dual-Simplex (Matlab R2015a).
% 
%  Usage: [v_e xm status]=p_exact_game(v,tol);
%
% Define variables:
%  output:
%  v_e      -- The exact game from v.
%  xm       -- The matrix of optimal imputation vectors. 
%  status   -- Returns 1 if the optimization terminated successfully. 
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
%   05/27/2024        1.9.2           hme
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
A1(N+1,:)=PlyMat(end,:);
A2=sparse(A1);
B1=[-v';v(N)];
ub=inf(1,n);
lb=vi';
v_e=zeros(1,N);
xm=zeros(N,n);
%opts=optimset('Simplex','on','LargeScale','off','MaxIter',256);
opts.Display='off';
opts.Simplex='on';
%opts.ActiveSet='on';
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
%% for dual-simplex
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);

parfor S=1:N
    C=PlyMat(S,:);
    try
       [xmin,fmin,exitflag]=linprog(C,A2,B1,[],[],lb,ub,opts);
    catch
       [xmin,fmin,exitflag]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
    end
    xm(S,:)=xmin';
    v_e(S)=fmin;
    if exitflag ~= 1 
       warning('Solution not optimal!')
    end
end    



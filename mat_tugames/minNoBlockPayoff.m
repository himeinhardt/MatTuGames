function MNBP=minNoBlockPayoff(v,tol);
% EXACT_GAME computes the minimum no blocking payoff from game v using 
% the Matlab's Optimization toolbox. Uses Dual-Simplex (Matlab R2015a).
% A core is non-empty iff mnbp<=v(N) 
%
%  Source: J. Zhao (2001), The relative interior of base polyhedron and the core, Economic Theory 18, 635-648. 
% 
%
%  Usage: NMBP=minNoBlockPayoff(v,tol);
%
% Define variables: 
%  Output of structure element NMBP:
%  mnbp     -- The minimum no blocking payoff.
%  pm       -- The minimum payoff allocation. 
%  crQ      -- Core is non-empty if 1 otherwise empty.
%  status   -- Returns 180 if the optimization terminated successfully. 
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
%   05/12/2016        0.8             hme
%                
if nargin<2
 tol=10^6*eps;
end

N=length(v);
[~,n]=log2(N);


k=1:n;
vi=v(bitset(0,k));

S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=-PlyMat;
A1(N,:)=[];
A2=sparse(A1);
B1=-v';
B1(N)=[];
ub=inf(1,n);
lb=vi';
pm=zeros(1,n);
%opts=optimset('Simplex','on','LargeScale','off','MaxIter',256);
opts.Display='off';
opts.Simplex='on';
%opts.ActiveSet='on';
opts.LargeScale='on';
opts.Algorithm='dual-simplex';
opts.TolFun=1e-10;
opts.TolX=1e-10;
opts.TolRLPFun=1e-10;
%% for dual-simplex
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);

C=ones(n,1);
[xmin,fmin,exitflag]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
pm=xmin';
mnbp=fmin;
if exitflag ~= 1 
   warning('Solution not optimal!')
end

crQ=mnbp<=v(N)+tol;

MNBP=struct('mnbp',mnbp,'pm',pm,'crQ',crQ,'status',exitflag);


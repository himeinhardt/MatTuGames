function [x1, fmin]=cplex_AntiPreNucl_llp(v,tol)
% CPLEX_ANTIPRENUCL computes the anti pre-nucleolus of game v using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
% 
%
% Usage: [x, alp]=cplex_AntiPreNucl_llp(v,tol)
% Define variables:
%  output:
%  x1        -- The anti pre-nucleolus of game v.
%  fmin      -- The maxmin excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/24/2017        0.9             hme
%   04/04/2020        1.9             hme
%                



if nargin<2
 tol=10^8*eps;
end
%tol=-tol;

N=length(v);
[~, n]=log2(N);
% solver parameter
ra = reasonable_outcome(v);
ub=[ra,inf];
x1=[];
lb=[-inf(1,n),-inf];
x0=[];
warning('off','all');
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else
%  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
  options.largescale='on';
  options.algorithm='dual-simplex';
  options.tolfun=1e-10;
  options.tolx=1e-10;
  options.tolrlpfun=1e-10;
  %%%% for dual-simplex
  % opts.MaxTime=9000;
  options.preprocess='none';
  options.tolcon=1e-6;
  options.maxiter=128;
  options.display='off';
  options.threads=3;
end
%warning('on','all');

S=1:N;
for k=1:n, A1(:,k) = bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[v';-v(N)];
C=[zeros(n,1);-1];

bA=find(A1(:,end)==0)';
while 1
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(lambda.ineqlin>tol))';
  bS1(end)=[];
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     warning('on','all');
     break;
  end
  it=0:-1:1-n;
  bA=[bA,bS2];
  mS2=rem(floor(bA(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+fmin;
  if rk==n 
     x=(-mS2\B1(bA))';
     warning('on','all');
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
end

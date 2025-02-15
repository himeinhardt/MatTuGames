function [x1, fmin]=cplex_cs_nucl(clv,cs,tol)
% CPLEX_CS_NUCL computes the nucleolus of game v w.r.t coalition structure cs 
% using cplexmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
% 
%
% Usage: [x, fmin]=clv.cplex_cs_nucl(cs,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game v with cs.
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  cs        -- A coalition structure provided as partition of N like [1 6].
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/23/2017        0.9             hme
%   04/04/2020        1.9             hme
%                


warning('off','all');
if nargin<3
 tol=10^8*eps;
end
%tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
vi=clv.tuvi;

if N==3
  x1=clv.StandardSolution();
  return
end


% solver parameter
ra = clv.reasonable_outcome();
k=1:n;
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if sum(vi)>v(N)
   error('sum of lower bound exceeds value of grand coalition! No solution can be found that satisfies the constraints.')
end
lb=[vi,-Inf]';
ub=[ra,Inf]';

x0=[];
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
  options.maxiter=10*(N+n);
  options.display='off';
  options.threads=3;
end
if iscell(cs)
   cs=clToMatlab(cs);
else
  cs=double(cs);
end
lcs=length(cs);

S=1:N;
for k=1:n, A1(:,k) = -bitget(S,k);end
a1=A1(cs,:);
A1(N+1:N+lcs,:)=-a1;
A1(:,n+1)=-1;
A1(N,:)=[];
A1(N:N+lcs-1,end)=0;
A2=sparse(A1);
vc=v;
vc(N)=[];
Bv=v(cs)';
B1=[-vc';Bv];
C=[zeros(n,1);1];

while 1
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(lambda.ineqlin>tol))';
  bS1(end)=[];
  bA=find(A1(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     warning('on','all');
     break;
  end
  it=0:-1:1-n;
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>-tol);
  if exitflag ~=1
     warning('on','all');
     warning('Prn:Exit','Probably no nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     warning('on','all');
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
end

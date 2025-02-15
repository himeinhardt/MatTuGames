function [x1, fmin]=glpk_cs_nucl(v,cs,tol)
% GLPK_CS_NUCL computes the nucleolus of game v w.r.t. coalition structure cs
% using glpkmex.
% 
% http://www.gnu.org/software/glpk/glpk.html
%
%
% Usage: [x, alp]=glpk_cs_nucl(v,cs,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game v.
%  fmin      -- The minmax excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  cs       -- A coalition structure provided as partition of N like [1 6].
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/31/2017        0.9             hme
%                



if nargin<3
 tol=10^6*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);

% solver parameter
ra = reasonable_outcome(v);
k=1:n;
vi=v(bitset(0,k));
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
if sum(vi)>v(N)
   error('sum of lower bound exceeds value of grand coalition! No solution can be found that satisfies the constraints.')
end

lb=[vi,-Inf];
ub=[ra,Inf];
%lb=[];
ctype=[];
vartype=[];
s=1; % minimization problem 
param.lpsolver=1; % simplex method

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
C=[zeros(1,n),1];


while 1
  [xmin, fmin, status, extra] = glpk(C,A1,B1,lb,ub,ctype,vartype,s,param);
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(extra.lambda<tol))';
  bS1(end)=[];
  bA=find(A1(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  it=0:-1:1-n;
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>tol);
  if status ~=5
     warning('Prn:Exit','Probably no nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
end

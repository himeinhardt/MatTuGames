function [x1, fmin]=glpk_nucl(clv,tol)
% GLPK_NUCL computes the nucleolus of game v using glpkmex.
% 
% http://www.gnu.org/software/glpk/glpk.html
%
%
% Usage: [x, alp]=glpk_nucl(clv,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game clv.
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/21/2012        0.3             hme
%                



if nargin<2
 tol=10^6*eps;
end
tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
vi=clv.tuvi;
if essQ==0
   error('Game is not essential!');
end
if N==3
  x1=clv.StandardSolution();
  fmin=-inf;
  return
end

% solver parameter
ra = clv.reasonable_outcome();
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end

lb=[vi,-Inf];
ub=[ra,Inf];
ctype=[];
vartype=[];
s=1; % minimization problem 
param.lpsolver=1; % simplex method


S=1:N;
A1=zeros(N,n);
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[-v';v(N)];
C=[zeros(1,n),1];
it=0:-1:1-n;
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

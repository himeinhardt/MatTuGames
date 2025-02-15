function [x1, fmin]=glpk_AntiPreNucl_llp(clv,tol)
% GLPK_ANTIPRENUCL computes the anti pre-nucleolus of game v using glpkmex.
% 
% http://www.gnu.org/software/glpk/glpk.html
%
%
% Usage: [x, alp]=clv.glpk_AntiPreNucl_llp(tol)
% Define variables:
%  output:
%  x1        -- The anti pre-nucleolus of game v.
%  fmin      -- The maxmin excess value.
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
%   08/05/2016        0.9             hme
%   03/25/2021        1.9             hme
%                



if nargin<2
 tol=10^6*eps;
end
tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
if N==3
  x1=clv.StandardSolution();
  return
end

% solver parameter
ra = clv.reasonable_outcome;
ub=[ra,Inf];
lb=[-inf(1,n),-Inf];
%lb=[];
ctype=[];
vartype=[];
s=1; % minimization problem 
param.lpsolver=1; % simplex method


S=1:N;
for k=1:n, A1(:,k) = bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=1;
A1(N:N+1,end)=0;
A2=sparse(A1);
B1=[v';-v(N)];
C=[zeros(1,n),-1];
bA=find(A1(:,end)==0)';
while 1
  [xmin, fmin, status, extra] = glpk(C,A2,B1,lb,ub,ctype,vartype,s,param);
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(extra.lambda<tol))';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  it=0:-1:1-n;
  bA=[bA,bS2];
  mS2=rem(floor(bA(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+fmin;
  if rk==n 
     x=(-mS2\B1(bA))';
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
end

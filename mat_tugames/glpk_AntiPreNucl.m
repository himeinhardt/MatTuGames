function [x1, fmin]=glpk_AntiPreNucl(v,tol)
% GLPK_ANTIPRENUCL computes the anti pre-nucleolus of game v using glpkmex.
% 
% http://www.gnu.org/software/glpk/glpk.html
%
%
% Usage: [x, alp]=glpk_AntiPreNucl(v,tol)
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
%   08/29/2014        0.5             hme
%                



if nargin<2
 tol=10^6*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);

% solver parameter
ra = reasonable_outcome(v);
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

while 1
  [xmin, fmin, status, extra] = glpk(C,A2,B1,lb,ub,ctype,vartype,s,param);
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
     warning('Prn:Exit','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
end

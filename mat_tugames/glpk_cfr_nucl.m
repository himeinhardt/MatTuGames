function [x1, fmin]=glpk_cfr_nucl(v,F,tol)
% GLPK_CFR_NUCL computes the nucleolus of game v with coalition formation restrictions 
% using glpkmex.
% 
% http://www.gnu.org/software/glpk/glpk.html
%
%
% Usage: [x, alp]=glpk_cfr_nucl(v,F,tol)
% Define variables:
%  output:
%  x1        -- The nucleolus of game vF.
%  fmin      -- The minmax excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  F        -- For instance, a characterization set for the nucleolus.
%              F must contain the grand coalition N.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/25/2017        0.9             hme
%                



if nargin<2
   error('A collection of sets F is required!');
elseif nargin<3
 tol=10^8*eps;
end

tol=-tol;

N=length(v);
[~, n]=log2(N);
k=1:n;
si=bitset(0,k);
%% Addining the singleton coalitions to F, 
%% since v is definded over F ans si.
F=unique([F,si]);
lf=length(F);

S=1:N;
lfNq=F(end)~=N;
if lfNq
   F(end+1)=N;
   lf=lf+1;
end
CS=S(ismember(S,F)==0);
vF=v;
vF(CS)=[];

if lfNq
   vF(end)=0;
end


% solver parameter
ra = reasonable_outcome(v);
vi=v(si);
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


for k=1:n, A1(:,k) = -bitget(F,k);end
A1(lf+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(lf:lf+1,end)=0;
A2=sparse(A1);
B1=[-vF';vF(lf)];
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
  bA(end)=[];
  bA=F(bA);
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

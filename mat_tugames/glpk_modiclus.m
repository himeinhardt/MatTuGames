function [x1, fmin]=glpk_modiclus(v,tol)
% GLPK_MODICLUS computes the modiclus of game v using glpkmex.
% 
% http://www.gnu.org/software/glpk/glpk.html
%
%
% Usage: [x, alp]=glpk_modiclus(v,tol)
% Define variables:
%  output:
%  x1        -- The modiclus of game v.
%  fmin      -- The minmax bi-excess value.
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
%   12/19/2017        0.9             hme
%                



if nargin<2
 tol=10^8*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);
dv=dual_game(v);

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
N1=N+1;
n1=2*n;
N2=2^n1-1;
geQ=all(v<=dv);
for k=1:n, A1(:,k) = -bitget(S,k);end
ii=1;
vs=zeros(N2,1);
A2=zeros(N2,n);
for k=1:N1;
    for jj=1:N1
        if geQ
          if k>1 && jj >1
             ii=(k-1)+(jj-1)*N1;
             A2(ii,:)=A1(k-1,:)+A1(jj-1,:);
             vs(ii)=v(k-1)+dv(jj-1);
           elseif k==1 && jj >1
             ii=N1*(jj-1);
             A2(ii,:)=A1(jj-1,:);
             vs(ii)=dv(jj-1);
           elseif k>1 && jj==1
             ii=k-1;
             A2(ii,:)=A1(k-1,:);
             vs(ii)=v(k-1);
           end
        else
          if k>1 && jj >1
             ii=N1*(k-1)+(jj-1);
             A2(ii,:)=A1(k-1,:)+A1(jj-1,:);
             vs(ii)=v(k-1)+dv(jj-1);
           elseif k==1 && jj >1
             ii=jj-1;
             A2(ii,:)=A1(jj-1,:);
             vs(ii)=dv(jj-1);
           elseif k>1 && jj==1
             ii=N1*(k-1);
             A2(ii,:)=A1(k-1,:);
             vs(ii)=v(k-1);
           end
        end
    end;
end


A2(N2+1,:)=-A2(end,:);
A2(:,end+1)=-1;
A2(N2:N2+1,end)=0;
B1=[-vs;v(N)+dv(N)];

%A2=sparse(A2);
C=[zeros(1,n),1];
it=0:-1:1-n1;
bA=find(A2(:,end)==0)';
while 1
  A2=sparse(A2);
  [xmin, fmin, status, extra] = glpk(C,A2,B1,lb,ub,ctype,vartype,s,param);
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(extra.lambda<tol))';
  bS1(end)=[];
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  bA=[bA,bS2];
  mS2=rem(floor(bA(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+fmin;
  if rk==n1
     x=(-mS2\B1(bA))';
     break;
  end
  A2(bS2,end)=0;
end

function [x1, fmin]=cplex_AntiNucl_llp(clv,tol)
% CPLEX_ANTINUCL_LLP computes the anti nucleolus of game v using glpkmex.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.8.0 and higher)
% 
%
% Usage: [x, fmin]=clv.cplex_AntiNucl_llp(tol)
% Define variables:
%  output:
%  x1        -- The anti nucleolus of game v.
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
%   02/09/2017        0.9             hme
%   02/24/2018        0.9             hme
%                


if nargin<2
 tol=10^8*eps;
end
%tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
vi=clv.tuvi;
if essQ==1
   error('Game is not anti essential!')
end


% solver parameter
ra = clv.smallest_amount();
k=1:n;
vi=v(bitset(0,k));
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
lb=[ra,-Inf]';
ub=[vi,Inf]';

x0=[];
warning('off','all');
%cplex = Cplex('null');
%str=cplex.getVersion;
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else 
  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
end
warning('on','all');
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

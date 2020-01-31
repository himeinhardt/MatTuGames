function [x1, alp]=CddModiclus(v,tol)
% CDDMODICLUS computes the modiclus of game v using cddmex.
% 
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% and the Matlab interface
% to the cdd solver (cddmex) http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
% Usage: [x, alp]=CddModiclus(v,tol)
% Define variables:
%  output:
%  x1        -- The modiclus of game v.
%  alp       -- The minmax bi-excess value.
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
%   12/21/2012        0.9             hme
%                



if nargin<2
 tol=10^6*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);
dv=dual_game(v);
ra = reasonable_outcome(v)';
vi = smallest_amount(v)';
S=1:N;
N1=N+1;
n1=2*n;
N2=2^n1-1;
geQ=all(v<=dv);
for k=1:n, A1(:,k) = -bitget(S,k);end
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
A2=[A2;-eye(n);eye(n)];
A2(:,end+1)=-1;
A2(N2:N2+1,end)=0;
B1=[-vs;v(N)+dv(N);vi;ra];
A1=A2;
objective=[zeros(1,n),1];
bA=find(A1(:,end)==0)';
it=0:-1:1-n1;
while 1
  IN=struct('obj',objective,'A',A1,'B',B1);
  OUT = cddmex('solve_lp_DS',IN);
  x=OUT.xopt;
  x1=x';
  x1(end)=[];
  alp=OUT.objlp;
  bS1=(find(OUT.lambda<tol))';
  bS1(end)=[];
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  bA=[bA,bS2];
  mS2=rem(floor(bA(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+alp;
  if rk==n1
     x=(-mS2\B1(bA))';
     break;
  end
  A1(bS2,end)=0;
end

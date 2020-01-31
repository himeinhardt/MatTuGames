function [x1, alp]=CddAntiNucl_llp(clv,tol)
% CDDANTINUCL_LLP computes the anti nucleolus of game v using cddmex.
% 
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% and the Matlab interface
% to the cdd solver (cddmex) http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Usage: [x, alp]=clv.CddAntiNucl_llp(tol)
% Define variables:
%  output:
%  x1        -- The anti nucleolus of game v.
%  alp       -- The minmax excess value.
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
%                



if nargin<2
 tol=10^6*eps;
end
tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
vi=clv.tuvi';
if essQ==1
   error('Game is not anti essential!')
end

ra = clv.smallest_amount()';
k=1:n;
vi=v(bitset(0,k))';
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end

S=1:N;
for k=1:n, A1(:,k) = bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1=[A1;-eye(n);eye(n)];
A1(:,end+1)=1;
A1(N:end,end)=0;
B1=[v';-v(N);ra;vi];
objective=[zeros(1,n),-1];
bA=find(A1(:,end)==0)';
while 1
  IN=struct('obj',objective,'A',A1,'B',B1);
  OUT = cddmex('solve_lp_DS',IN);
  x=OUT.xopt;
  x1=x';
  x1(end)=[];
  alp=OUT.objlp;
  bS1=(find(OUT.lambda(1:N+1)<tol))';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  it=0:-1:1-n;
  bA=[bA,bS2];
  mS2=rem(floor(bA(:)*pow2(it)),2);
  rk=rank(mS2);
  B1(bS2)=B1(bS2)+alp;
  if rk==n
     x=(-mS2\B1(bA))';
     break;
  end
  A1(bS2,end)=0;
end

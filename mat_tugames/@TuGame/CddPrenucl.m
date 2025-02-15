function [x1, alp]=CddPrenucl(clv,tol)
% CDDPRENUCL computes the pre-nucleolus of game clv using cddmex.
% 
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% and the Matlab interface
% to the cdd solver (cddmex) http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Usage: [x, alp]=CddPrenucl(clv,tol)
% Define variables:
%  output:
%  x1        -- The pre-nucleolus of game clv.
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
%   12/12/2012        0.3             hme
%                



if nargin<2
 tol=10^8*eps;
end
tol=-tol;

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
if N==3
  x1=clv.StandardSolution();
  alp=-inf;
  return
end
S=1:N;
A1=zeros(N,n);
for k=1:n, A1(:,k) = -bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(N:N+1,end)=0;
B1=[-v';v(N)];
objective=[zeros(1,n),1];
it=0:-1:1-n;
while 1
  IN=struct('obj',objective,'A',A1,'B',B1);
  OUT = cddmex('solve_lp_DS',IN);
  x=OUT.xopt;
  x1=x';
  x1(end)=[];
  alp=OUT.objlp;
  bS1=(find(OUT.lambda<tol))';
  bA=find(A1(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1,
     break;
  end
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>tol);
  if OUT.how ~=1
     warning('Prn:Exit','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  B1(bS2)=B1(bS2)+alp;
end

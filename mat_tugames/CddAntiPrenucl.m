function [x1, alp]=CddAntiPrenucl(v,tol)
% CDDANTIPRENUCL computes the anti pre-nucleolus of game v using cddmex.
% 
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% and the Matlab interface
% to the cdd solver (cddmex) http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Usage: [x, alp]=CddAntiPrenucl(v,tol)
% Define variables:
%  output:
%  x1        -- The anti pre-nucleolus of game v.
%  alp       -- The maxmin excess value.
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
 tol=10^8*eps;
end
tol=-tol;

N=length(v);
[~, n]=log2(N);
S=1:N;
for k=1:n, A1(:,k) = bitget(S,k);end
A1(N+1,:)=-A1(end,:);
A1(:,end+1)=1;
A1(N:N+1,end)=0;
B1=[v';-v(N)];
objective=[zeros(1,n),-1];

while 1
  IN=struct('obj',objective,'A',A1,'B',B1);
  OUT = cddmex('solve_lp_DS',IN);
  x=OUT.xopt;
  x1=x';
  x1(end)=[];
  alp=OUT.objlp;
  bS1=(find(OUT.lambda<tol))';
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
  if OUT.how ~=1
     warning('Prn:Exit','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     break;
  end
  A1(bS2,end)=0;
  B1(bS2)=B1(bS2)+alp;
end

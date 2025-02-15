function mwgs=MinimalRep(sv,th);
% MINIMALREP computes from a simple game v and a threshold th of the 
% minimal representation of an homogeneous weighted majority game.
%
% Source: Sudhoelter (1996)
%
% Usage: mwgs=NetworkMinimalRep(sv,th);
%
% Define variables:
%  output:
%  mwgs     -- Vector of minimal weights.
%
%  input:
%  sv       -- A simple game of length 2^n-1.
%  th       -- Threshold to pass a bill (positive number).
%  tol      -- Numerical tolerance.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/26/2017        0.9             hme
%                    

if nargin<3
   tol=10^8*eps;
end
%[mW,wC]=getMinimalWinning(sv);
mW=getMinimalWinning(sv);
N=length(sv);
[~, n]=log2(N);

for kk=1:n
    M(:,kk)=bitget(mW,kk)==1;
end
dM=double(M);
[~,S,~] = svd(dM);
s = diag(S);
ls=s>tol;
m2=nnz(ls);
M=M(ls,:);
M=M(:,ls);
b=ones(m2,1)*th;
x=pinv(double(M))*b;
lx=length(x);
if lx<n
   mwgs=x';
   mwgs(n)=0;
else
   mwgs=x';
end

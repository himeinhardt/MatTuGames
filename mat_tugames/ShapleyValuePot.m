function shv=ShapleyValuePot(v)
% SHAPLEY_VALUEPOT computes the Shapley-value of a TU-game v from 
% the relevant parts of the potential. 
%
% Usage: shv=ShapleyValuePot(v)
%
% Define variables:
%  output:
%  sh       -- The Shapley-value of a TU-game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/04/2014        0.5             hme
%                

N=length(v);
[~, n]=log2(N);
if N==1
  shv=v;return;
 else
end

int=1-n:1:0;
k=1:n;
shv=zeros(1,n);
Nk=bitset(N,k,0);
[sNk,idx]=sort([Nk N]);
lNk=n+1;
% Try this, if you have a memory problem.
% Might be slow!!
%nzm=3*N;
%pb=spalloc(N,lNk,nzm);
pb=zeros(N,lNk);
for k=1:lNk
   sS=Subsets(sNk(k),n);
   mat=rem(floor(sS(:)*pow2(int)),2);
   cS=mat*ones(n,1);
   ck=rem(floor(sNk(k)*pow2(int)),2);
   sk=ck*ones(n,1);
   pb(sS,k)=(factorial(cS-1).*factorial(sk-cS))/factorial(sk);
end
%pb=sparse(pb);
pot=v*pb(:,idx);
idx(end)=[];
shv(idx)=pot(end)-pot(idx);

%--------------------------------------
function sS=Subsets(S,n)

it=0:-1:1-n;
vecS=rem(floor(S(:)*pow2(it)),2);

J=1:n;
slcP=vecS==0;
sP=J(slcP);

S1=1:S; 

if (2^n-1)==S
  sS=S1;
else 
 lsP=length(sP);
 Tni=cell(lsP);
 for k=1:lsP
  Tni{k}=bitget(S1,sP(k))==0;
 end

 cls=size(Tni);
 ls1=length(S1);
 R=true(1,ls1);
 for k=1:cls(:,2)
  R=Tni{k} & R;
 end
 sS=S1(R);
end

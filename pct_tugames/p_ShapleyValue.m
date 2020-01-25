function shv=p_ShapleyValue(v)
% P_SHAPLEY_VALUE computes the Shapley-value of a TU-game v and 
% the relevant parts of the potential using Matlab's PCT.
%
% Usage: shv=p_ShapleyValue(v)
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
elseif iscolumn(v)
    v=v';
end


S=1:N;
parfor k=1:n
  mat(:,k)=bitget(S,k)==1;
end
clS=mat*ones(n,1);
clear mat S;
k=1:n;
Nk=N-2.^(k-1);
[sNk,idx]=sort([Nk N]);
lNk=n+1;
ck=zeros(1,lNk);
ck(1:n)=n-1;
ck(n+1)=n;
pot=zeros(1,lNk);
parfor k=1:lNk
   T=sNk(k);
   sS=Subsets(T,n);
   vS=v(sS);
   cS=clS(sS);
   sk=ck(k);
   pb=(factorial(cS-1).*factorial(sk-cS))/factorial(sk);
   pot(k)=vS*pb;
end
idx(lNk)=[];
shv=pot(lNk)-pot(idx);

%--------------------------------------
function sS=Subsets(S,n)

it=0:-1:1-n;
slcP=rem(floor(S(:)*pow2(it)),2)==0;

J=1:n;
sP=J(slcP);

S1=1:S;

if (2^n-1)==S
  sS=S1;
else
 lsP=length(sP);
 ls1=length(S1);
 Tni=false(lsP,ls1);
 for k=1:lsP
  Tni(k,:)=bitget(S1,sP(k))==0;
 end
 R=true(1,ls1);
 for k=1:lsP
  R=Tni(k,:) & R;
 end
 sS=S1(R);
end

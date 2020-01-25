function x=LS_Nucl(v)
% LS_NUCL computes the least square nucleolus of a game.
%
% Usage: x=LS_Nucl(v)
% Define variables:
%  output:
%  x        -- Least square nucleolus 
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/10/2013        0.4             hme
%     
    
N=length(v);
[~, n]=log2(N);

J=1:n;
vi=v(bitset(0,J));
slb=sum(vi)>v(N);
if slb==1
   error('Game is not essential!')
end

    
S=1:N;
a=cell(1,n);
sV=zeros(1,n);

for k=1:n
    a{k}=S(bitget(S,k)==1);
    sV(k)=sum(v(a{k}));
end     
nV=sum(sV);
sN=2^(n-2);
x=(v(N)/n)+(n*sV-nV)/(n*sN);

M=J(x<0);
if isempty(M)
   return; 
end    
while 1
  M0=M;  
  CM=setdiff(J,M);
  if isempty(CM)==0
     xm=sum(x(M));
     m=length(M);
     x(CM)=x(CM)+xm/(n-m);
     x(M)=0;
  end 
  M2=J(x<0);
  M=unique([M,M2]);
  if M==M0
     break; 
  end    
end  

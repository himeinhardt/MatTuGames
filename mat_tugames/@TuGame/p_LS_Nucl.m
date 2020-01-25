function x=p_LS_Nucl(clv)
% P_LS_NUCL computes the least square nucleolus of a game
% using Matlab's PCT.
%
% Usage: x=p_LS_Nucl(clv)
% Define variables:
%  output:
%  x        -- Least square nucleolus 
%
%  input:
%  clv      -- TuGame class object.
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
    
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
if essQ==0
   error('Game is not essential!')
end



S=1:N;
sV=zeros(1,n);

 parfor k=1:n
     a=S(bitget(S,k)==1);
     sV(k)=sum(v(a));
 end     
nV=sum(sV);
sN=2^(n-2);
x=(v(N)/n)+(n*sV-nV)/(n*sN);

J=1:n;
M=J(x<0);
if isempty(M)
   return; 
end    
while 1
  M0=M;  
  CM=setdiff(J,M);
  if isempty(CM)==0
     xm=sum(x([M]));
     m=length(M);
     x([CM])=x([CM])+xm/(n-m);
     x([M])=0;
  end 
  M2=J(x<0);
  M=unique([M,M2]);
  if M==M0
     break; 
  end    
end  

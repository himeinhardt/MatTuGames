function [pot,shv,pb]=p_potential(v)
% P_POTENTIAL computes the potential of a TU-game from a game basis.
% Very slow, for n>14 this function needs some time to complete. If
% you are not interested in the game basis use the function
% Potential to compute the potential more efficiently.
%
% Usage: [pot,shv,pb]=p_potential(v)
%
% Define variables:
%  output:
%  pot      -- The potential of the game v.
%  sh       -- The Shapley-value of a TU-game v.
%  pb       -- Game basis.
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
pot=zeros(1,N);
S=1:N;
parfor k=1:n
  mat(:,k)=bitget(S,k)==1;
end
clS=mat*ones(n,1);
clear mat S;


if nargout < 3
    parfor k=1:N
      pb=zeros(N,1);  
      sS=Subsets(k,n);
      cS=clS(sS);
      sk=clS(k);
      if sk==cS
         pb(sS)=1; 
      else    
         pb(sS)=(factorial(cS-1).*factorial(sk-cS))/factorial(sk);
      end
      pot(k)=v*pb;
    end
  k=1:n;
  Nk=N-2.^(k-1);
  shv=pot(N)-pot(Nk);    
else
  pb=zeros(N);
  parfor k=1:N
      mytmp = zeros(N,1);
      sS=Subsets(k,n);
      cS=clS(sS);
      sk=clS(k);
      if sk==cS
         mytmp(sS,:)=1;
      else
         mytmp(sS,:)=(factorial(cS-1).*factorial(sk-cS))/factorial(sk);
      end
      pot(k)=v*mytmp;
      pb(:,k)=mytmp;
  end    
  k=1:n;
  Nk=N-2.^(k-1);
  shv=pot(N)-pot(Nk);
end


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

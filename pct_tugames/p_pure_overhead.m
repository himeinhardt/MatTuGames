function pv=p_pure_overhead(n)
% P_PURE_OVERHEAD computes the pure overhead games.
% 
% Usage: pv=p_pure_overhead(n)
%

% Define variables:
%  output:
%  pv       -- Matrix of pure overhead games.
%
%  input:
%  n        -- number of players involved.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/03/2013        0.4             hme
%                


N=2^n-1;
lb=zeros(N);
S=1:N;
int=1-n:1:0;
pv=zeros(N);

parfor k=1:N
   aS=bitand(S,k);
   lS=length(aS);
   idx1=S(aS==0);
   idx2=S(aS>0);
   aS(idx1)=[];
   mat=rem(floor(aS(:)*pow2(int)),2)';
   temp = zeros(1,lS);
   temp(idx2)=ones(1,n)*mat-1;  
   pv(k,:)=temp;
end

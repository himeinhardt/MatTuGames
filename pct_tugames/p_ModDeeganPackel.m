function mdpidx=p_ModDeeganPackel(th,w_vec)
% P_ModDeeganPackel computes the modified Deegan-Packel index from the set of winning coalitions
% using Matlab's PCT.
%
% Usage: dpidx=p_DeeganPackel(th,w_vec)
% Define variables:
%  output:
%  mdpidx   -- The modified Deegan-Packel index.
%
%  input:
%  th       -- Threshold to pass a bill (positive number).
%  w_vec    -- Vector of weights.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/01/2020        1.9             hme
%


n=length(w_vec);
sWCk=zeros(1,n);
[mWC,wC]=p_minimal_winning(th,w_vec);
%wC=winning_coalitions(mWC,n);
m=length(wC);
mdpidx=zeros(1,n);
parfor k=1:n;
  if w_vec(k)~=0 
     mWCk=wC(bitget(wC,k)==1);
     sk=length(mWCk);
     A=zeros(sk,1);
     for jj=1:n, 
         A(:,jj) = bitget(mWCk,jj);
     end
     ssk=sum(1./(A*ones(n,1)));
     mdpidx(k)=ssk/m;
  end   
end


function dpidx=DeeganPackel(th,w_vec)
% DeeganPackel computes the Deegan-Packel index from the set of minimal winning coalitions.
%
% Usage: dpidx=DeeganPackel(th,w_vec)
% Define variables:
%  output:
%  dpidx    -- The Deegan-Packel index.
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
mWC=minimal_winning(th,w_vec);
m=length(mWC);
dpidx=zeros(1,n);
for k=1:n;
  mWCk=mWC(bitget(mWC,k)==1);
  sk=length(mWCk);
  A=zeros(sk,1);
  for jj=1:n, 
      A(:,jj) = bitget(mWCk,jj);
  end
  ssk=sum(1./(A*ones(n,1)));
  dpidx(k)=ssk/m; 
end


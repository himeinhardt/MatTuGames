function hidx=p_holler(th,w_vec)
% P_HOLLER computes the Holler index from the set of minimal winning coalitions
% while using Matlab's PCT.
%
% Usage: hidx=p_holler(th,w_vec)
%
% Define variables:
%  output:
%  hidx      -- The Holler index.
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
%   05/25/2014        0.5             hme
%

n=length(w_vec);
sWCk=zeros(1,n);
mWC=p_minimal_winning(th,w_vec);
parfor k=1:n;
  mWCk=mWC(bitget(mWC,k)==1);
  sWCk(k)=length(mWCk);
end
hidx=sWCk./(sWCk*ones(n,1));


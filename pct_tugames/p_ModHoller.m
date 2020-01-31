function mhidx=p_ModHoller(th,w_vec)
% P_MODHOLLER computes the modified Holler index from the set of winning coalitions
% while using Matlab's PCT.
%
% Usage: mhidx=p_ModHoller(th,w_vec)
% Define variables:
%  output:
%  mhidx     -- Modified Holler index.
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
%   09/25/2018        1.0             hme
%

n=length(w_vec);
sWCk=zeros(1,n);
mWC=p_minimal_winning(th,w_vec);
WC=winning_coalitions(mWC,n);
parfor k=1:n;
  mWCk=WC(bitget(WC,k)==1);
  sWCk(k)=length(mWCk);
end
mhidx=sWCk./(sWCk*ones(n,1));


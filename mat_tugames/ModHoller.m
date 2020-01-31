function mhidx=ModHoller(th,w_vec)
% MODHOLLER computes a modified  Holler index from the set of winning coalitions.
% This avoids the violation of local monotonicity (Holler, 2018).
%
% Usage: mhidx=ModHoller(th,w_vec)
% Define variables:
%  output:
%  hidx      -- Modified Holler index.
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
%   09/22/2018        1.0             hme
%


n=length(w_vec);
sWCk=zeros(1,n);
mWC=minimal_winning(th,w_vec);
WC=winning_coalitions(mWC,n);
for k=1:n;
  mWCk=WC(bitget(WC,k)==1);
  sWCk(k)=length(mWCk);
end
mhidx=sWCk./(sWCk*ones(n,1));


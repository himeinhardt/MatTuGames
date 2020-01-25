function [mW wC wm]=p_minimal_winning(th,w_vec)
% P_MINIMAL_WINNING computes from the threshold th and
% the weights w_vec the minimal winning coalitions 
% while using Matlab's PCT.
%
% Usage: [mW wC wm]=p_minimal_winning(th,w_vec);
%
% Define variables:
%  output:
%  mW       -- The list/vector of minimal winning coalitions.
%  wC       -- The list of winning coalitions. 
%  wm       -- Weighted majority game of length 2^n-1.
%  input:
%  th       -- Threshold to pass a bill (positive number).
%  w_vec    -- Vector of weights.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/29/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%                    

wm=weighted_majority(th,w_vec);
N=length(wm);
[~, n]=log2(N);
S=1:N;
wC=S(wm);
lwC=length(wC);

parfor k=1:lwC
  sS=SubSets(wC(k),n);
  [vl, idx]=max(wm(sS));
  smS=sS(idx);
  if smS==wC(k)
      mWS(k)=wC(k);
  end
end

mW=mWS(logical(mWS>0));

function mpgi=ModPGI(th,w_vec)
% MODGPI computes a modified  public good index from the set of winning coalitions.
% This avoids the violation of local monotonicity (Holler, 2018).
%
% Usage: mpgi=ModPGI(th,w_vec)
% Define variables:
%  output:
%  hidx      -- Modified Public Good/Holler index.
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
%   10/03/2020        1.9             hme
%

mpgi=ModHoller(th,w_vec);


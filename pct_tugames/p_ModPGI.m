function mpgi=p_ModPGI(th,w_vec)
% P_MODPGI computes the modified public good index from the set of winning coalitions
% while using Matlab's PCT. Same as modified Holler index.
%
% Usage: mpgi=p_ModPGI(th,w_vec)
% Define variables:
%  output:
%  mpgi     -- Modified Public Good/Holler index.
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

mpgi=p_ModHoller(th,w_vec);


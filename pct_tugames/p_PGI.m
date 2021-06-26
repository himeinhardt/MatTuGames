function pgi=p_GPI(th,w_vec)
% P_PGI computes the public good index from the set of minimal winning coalitions
% while using Matlab's PCT. Same as Holler index.
%
% Usage: pgi=p_PGI(th,w_vec)
%
% Define variables:
%  output:
%  pgi      -- The Public Good/Holler index.
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

pgi=p_holler(th,w_vec);


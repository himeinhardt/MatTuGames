function pgi=PGI(th,w_vec)
% PGI computes the public good index from the set of minimal winning coalitions.
% Same as Holler index (holler).
%
% Usage: gpi=PGI(th,w_vec)
% Define variables:
%  output:
%  pgi      -- The Holler/Public good index.
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


pgi=holler(th,w_vec);



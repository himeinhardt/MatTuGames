function v=weighted_majority(th,w_vec)
% WEIGHTED_MAJORITY computes from the threshold th and
% the weights w_vec a weighted majority game.
%
% Usage: v=weighted_majority(th,w_vec)
% Define variables:
%  output:
%  v        -- A weighted majority game.
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
%   08/06/2010        0.1 beta        hme
%   07/02/2012        0.2 beta        hme
%   09/21/2012        0.3 beta        hme
%                



n=length(w_vec);
sumw=w_vec(1); for ii=2:n, sumw=[sumw w_vec(ii) sumw+w_vec(ii)]; end
v=sumw>=th;


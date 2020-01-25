function v=greedy_bankruptcy(E,d_vec)
% GREEDY_BANKRUPTCY computes for a bankruptcy situation (E,d_vec)
% the corresponding greedy bankruptcy game.
%
% Usage: v=greedy_bankruptcy(E,d_vec)
% Define variables:
%  output:
%  v        -- A greedy TU-bankruptcy game. The dual of the modest
%              bankruptcy game.
%  input:
%  E        -- Estate E (positive number).
%  d        -- Vector of claims of the claimants.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/10/2010        0.1 beta        hme
%   07/02/2012        0.2 beta        hme
%   09/18/2012        0.2             hme
%                

n=length(d_vec);
sumd=d_vec(1); for k=2:n sumd=[sumd d_vec(k) sumd+d_vec(k)]; end
v=min(E,sumd);

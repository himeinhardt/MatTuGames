function v=bankruptcy_game(E,d_vec)
% BANKRUPTCY_GAME computes for a bankruptcy situation (E,d_vec)
% the corresponding bankruptcy game.
%
% Usage: v=bankruptcy_game2(E,d_vec)
% Define variables:
%  output:
%  v        -- A TU-bankruptcy game.
%  input:
%  E        -- Estate E (positive number s.t. E <= sum(d)).
%  d        -- Vector of claims of the claimants.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/06/2010        0.1 beta        hme
%   07/02/2012        0.2 beta        hme
%   09/17/2012        0.3 beta        hme
%   07/04/2014        0.4             hme
%   07/24/2018        1.0             hme
%                


n=length(d_vec);
alp=E-sum(d_vec);
Dm=d_vec(1); for k=2:n Dm=[Dm d_vec(k) Dm+d_vec(k)]; end
dE=alp+Dm;
v=max(dE,0);

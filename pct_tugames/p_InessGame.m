function vi=p_InessGame(x)
% P_INESSRGAME recovers from the Shapley-value/payoff of a default game the
% corresponding inessential game using Matlab's PCT. The inessential game vi has x as
% its Shapley value. Reverse operation of the function DecomposeGame(). 
% For n>12 this function needs some time to complete.
%
% Usage: vi=p_InessGame(x)
%
% Define variables:
%  output:     structure elements
%  vi       -- An inessential TU-game.
%
%  input:
%  x        -- A payoff vector/Shapley value of a default game.
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/15/2014        0.5             hme
%

    
n=length(x);
lb=p_linear_basis(n,'sparse');    
k=1:n;
sC=2.^(k-1);
Z=lb(:,sC);
vi=x*Z';

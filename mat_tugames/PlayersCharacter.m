function PC=PlayersCharacter(th,w_vec)
% PlayersCharacter partitions the set of players of the weighted majority game into
% the character Sum, Step, and Null-Player.
%
% Source:  J. Rosenmueller, Homogeneous Games: Recursive Structure and Computation,  
%          Mathematics of Operations Research, Vol. 12, No. 2 (May, 1987), pp. 309-33
%
%
% Usage: PC=PlayersCharacter(th,w_vec);
%
% Define structure variables:
%  output:
%  sums     -- Returns all players with character sum.
%  steps    -- Returns all players with character step.
%  nlp      -- Returns all players with character null-player.
%  th       -- Original threshold to pass a bill (positive number).
%  w_vec    -- Original vector of weights (descend ordering).
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
%   10/01/2021        1.9             hme
%          
PPly=PartitionPlySet(th,w_vec);
%% Partition of the Player set.
PC.sums=PPly.sums; % Sums 
PC.steps=PPly.steps; % Steps 
PC.npl=PPly.npl; % Null Player
PC.th=th;
PC.w_vec=w_vec;

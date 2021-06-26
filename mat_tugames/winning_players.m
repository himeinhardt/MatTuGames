function wpl=winning_players(th,w_vec)
% WINNING_PLAYERS computes from a pre-defined set of winning 
% coalitions (e.g. minimal winning coalitions) the set of winning players,
% i.e., player i for which v(i) = 1 holds.     
%
% Source: N. Megiddo (1994) editor, Essays in Game Theory: In Honor of M. Maschler, Chap. 13.
%    
% Usage: wpl=winning_players(th,w_vec)
%
% Define variables:
%  output:
%  wpl      -- The list of winning players. 
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
%   01/22/2021        01.9            hme
%                    

n=length(w_vec);
pl=1:n;
wpl=pl(w_vec(pl)>=th);


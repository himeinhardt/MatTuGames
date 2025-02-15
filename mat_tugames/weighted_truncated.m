function [wmg_tr,MR]=weighted_truncated(th,w_vec);
% WEIGHTED_TRUNCATED computes from the threshold th and
% the weights w_vec a truncated weighted majority game. 
% The game is obtained from the last player that is equivalent 
% to the first step. The game contains only steps.
%
% Source: Sudhoelter (1996), Star-shapedness of the kernel for homogeneous games.
%
% Usage: v=weighted_truncated(th,w_vec)
% Define variables:
%  output:
%  wmg_tr   -- A truncated weighted majority game (reduced game).
%
% Field variables of structure MR:
%  output:
%  mrQ      -- Returns one (true) whenever a minimal homogeneous representation was found, 
%              otherwise zero (false).
%  hrQ      -- Returns one whenever the representation is homogeneous.
%  mwgs     -- Vector of minimal weights.
%  th       -- Quorum of the minimal homogeneous representation.
%  eQ       -- Returns one whenever the both weighed majority games 
%              i.e., (th0,w_vec0) vs. (th,mwgs), are equal.
%  smpl     -- Players with character sum.
%  stpl     -- Players with character step.
%  npl      -- Players with character null-player.
%  th0      -- Original quorum.
%  wvec0    -- Original weights.
%
%  input:
%  th       -- Threshold/quorum to pass a bill (positive number).
%  w_vec    -- Vector of weights (descend ordering).
%  hrQ      -- Set one (true) for a homogeneous representation, 
%              otherwise zero (false) or empty set []. Suppresses 
%              the costly check of homogeneity. 
%  tol      -- Numerical tolerance.
%
%  input:
%  th       -- Threshold to pass a bill (positive number).
%  w_vec    -- Vector of weights (descend order).


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/22/2021        1.9             hme
%   06/23/2022        1.9.1           hme
%

w_vec=sort(w_vec,'descend');
PC=PlayersCharacter(th,w_vec);
stpl=PC.steps;
bd=stpl(end);
eqnpl=PlyEqFirstStep(th,w_vec);
if eqnpl.lp < bd  
   lm=th-sum(w_vec(eqnpl.lp+1:bd));
   ms=w_vec(1:eqnpl.lp);
   wmg_tr=weighted_majority(lm,ms);
   if nargout==2
      MR=min_homogrep(lm,ms,1); %% must be homogeneous!! Is assumed!!
   end
else
   msg01='No equivalent step found!';	
   warning('PEqF:Wr',msg01);	
   wmg_tr=[];
   if nargout==2
     MR=[];  
   end       
end	

function [hrQ hr wmg]=homogeneous_representationQ(th,w_vec)
% HOMOGENEOUS_REPRESENTATIONQ  checks whether the weighted majority game
% derived from (th,w_vec) possesses a homogeneous representation.
% Th. (Sudhoelter 1996): The pre-kernel of every homogeneous weighted
% majority game is star-shaped. 
%
%
% Usage: [hrQ hr wmg]=homogeneous_representationQ(th,w_vec);
%
% Define variables:
%  output:
%  hrQ      -- Returns 1 if the weighted majority game wmg has a
%              homogeneous representation, otherwise 0.
%  hr       -- The list/vector of minimal winning coalitions.
%  wmg      -- The underlying weighted majority game.
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
%   05/29/2012        0.2 beta        hme
%   10/17/2012        0.3             hme
%                    

[mW, ~, wmg]=minimal_winning(th,w_vec);
n=length(w_vec);
for k=1:n, Clm(:,k)=bitget(mW,k); end
qtvec=Clm*w_vec';
qS=qtvec==th;
hrQ=all(qS);
if hrQ
   hr=mW;
else
   hr=0;
end

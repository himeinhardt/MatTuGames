function sh=ShapAirPortMod2(C,nj)
% SHAPAIRPORTMOD2 computes the Shapley value from a airport capital cost problem.
% 
% Source: S. C. Littlechild and G. Owen. A simple expression for the shapely value in a special case. Management Science, 20(3):370–372, 1973. ISSN 00251909, 15265501. URL http://www.jstor.org/stable/2629727.
%    
% Usage: sh=ShapAirPortMod2(C,nj)
% Define variables:
%  output:
%  sh        -- The Shapley value of the associated airport capital cost game.
%
%  input:
%  C        -- Vector of capital costs per types of length m.
%  nj       -- Vector of plane movements per types of length m.
%  

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/29/2023        1.9.1             hme
%
msg=nargchk(2,2,nargin); % check for legal number of input arguments.
error(msg);

m=length(C); % Determining the length of the cost vector.
smNj=sum(nj);  % The overall plane movements.
cumMv=smNj-cumsum(nj); % Determining the cumulative plane movments of types. 
cumMv(m)=[]; % Deleting zero from the list.
cumMv=[smNj,cumMv]; % Getting the reversed order of the cumulative plane movements. 
tC2=C; % Copying the cost vector to a new variable.
tC2(m)=[]; % Deleting the cost of type m.
tC2=[0,tC2]; % Shifting the vector one step to the east.
dC=C-tC2; % Determining the cost increment for all types. 
y=dC./cumMv; % Dividing the cost increment of types by the cumulative plane movements.
sh=cumsum(y); % Determining the cumulative sum of the ratio for each type, this is the Shapley value. 
end % Resevered word end to indicate the closing of the M-function.

    

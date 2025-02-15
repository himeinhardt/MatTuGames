function sh=ShapleyAirportProb(C,nj)
% SHAPLEYAIRPORTPROB computes the Shapley value from a airport capital cost problem.
% 
% Source: S. C. Littlechild and G. Owen. A simple expression for the shapely value in a special case. Management Science, 20(3):370–372, 1973. ISSN 00251909, 15265501. URL http://www.jstor.org/stable/2629727.
%    
% Usage: sh=ShapleyAirportProb(C,nj)
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
    
m=length(C);
sh=zeros(1,m);
% Determining recursively the Shapley value.
for k=1:m 
    if k-1<1
       sh(k)=C(k)/sum(nj(k:m));
    else   
       sh(k)=sh(k-1)+(C(k)-C(k-1))/sum(nj(k:m));
    end   
end    

    
    

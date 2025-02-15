function tv=TauValAirportProb(C,nj)
% TAUVALAIRPORTPROB computes the tau value from a airport capital cost problem.
% 
% Source: Theo Driessen (1988), Cooperative Games, Solutions and Applications
%    
% Usage: tv=TauValAirportProb(C,nj)
% Define variables:
%  output:
%  tv       -- The tau value of the associated airport capital cost game.
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

% Pre-allocating output variables.            
m=length(C);
pl=1:m;
tv=zeros(1,m);
% Determining the tau value.
if nj(m)>=2
   tv=(C(m).*C)./sum(C.*nj);    
else
   C1=C;
   nj1=nj;
   C1(m)=[];
   nj1(m)=[];
   tv(1:m-1)=(C1(end).*C1)./(sum(C1.*nj1)+C1(end));
   tv(m)=tv(m-1)+(C(m)-C(m-1));
end

    

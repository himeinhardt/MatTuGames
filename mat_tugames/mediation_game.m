function v=mediation_game(n,k)
% MEDIATION_GAME(n,k) assigns its outside option to a coalition
% without player k (employer) the other players are the employees.
% Otherwise, the coalition receives its quadratic production capacity.
% Usage: v=mediation_game(n,k)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  n        -- An integer to specify the number of persons involved in
%              the game.
%  k        -- An integer to specify the employer.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/20/2013        0.4             hme
%                
N=2^n-1;
S=1:N;
csz=zeros(1,N);
v=zeros(1,N);
for ii=1:n, PlyMat(:,ii) = bitget(S,ii)==1;end
csz=PlyMat*ones(n,1);
lg1=PlyMat(:,k)';
lg2=lg1==0;
cs=S(lg1);
es=S(lg2);
v(cs)=(csz(lg1)-1).^2;
v(es)=csz(lg2);




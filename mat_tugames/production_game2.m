function v=production_game2(n,ct)
% production_GAME2(n) zeros, otherwise a constant minus a cardinality.
%
% Usage: v=production_game2(n,k)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%
% input: 
%  n        -- An integer to specify the number of persons involved in
%              the game.
%  ct       -- A constant, which must be an integer.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/01/2016        0.8             hme
%

if nargin<2
   ct=10;
%   sr=4;
end


N=2^n-1;
S=1:N;
csz=zeros(1,N);
v=zeros(1,N);
k=1:n;
for ii=1:n, PlyMat(:,ii) = bitget(S,ii) == 1;end

for kk=1:N 
    pk=k(PlyMat(kk,:));
    if length(pk)==1
       v(kk)=pk;
    elseif length(pk)==2
        if abs(pk(1)-pk(2))>0
           v(kk)=ct-pk(1)-pk(2);
        end
    elseif length(pk)==3
        v(kk)=ct;
    else 
        v(kk)=ct+length(pk);
    end
end
%v(N)=ct+sr*length;




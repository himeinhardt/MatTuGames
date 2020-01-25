function gb=game_basis(n)
% GAMES_BASIS computes the game basis for n-players.
% For n>14 this function needs some time to complete.
%
% Usage: gb=game_basis(n)
% Define variables:
%  output:
%  gb       -- Game basis.
%
%  input:
%  n        -- Number of players.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/09/2010        0.1 beta        hme
%   10/27/2012        0.3             hme
%   05/15/2014        0.5             hme
%                

if nargin<1
  error('At least the number of player must be given! This must be a positive number.')
 else
end
if length(n)>1
  error('The input argument n must be an integer!')
 elseif n<=0
  error('The input argument n must be a positive number!')
 else
end

N=2^n-1;
gb=false(N);
S=1:N;

for k=1:N
   sS=Subsets(k,n);
   gb(k,:)= ismembc(S,sS);
end



%--------------------------------------
function sS=Subsets(S,n)

it=0:-1:1-n;
slcP=rem(floor(S(:)*pow2(it)),2)==0;

J=1:n;
sP=J(slcP);

S1=1:S; 

if (2^n-1)==S
  sS=S1;
else 
 lsP=length(sP);
 Tni=cell(lsP);
 for k=1:lsP
  Tni{k}=bitget(S1,sP(k))==0;
 end

 cls=size(Tni);
 ls1=length(S1);
 R=true(1,ls1);
 for k=1:cls(:,2)
  R=Tni{k} & R;
 end
 sS=S1(R);
end

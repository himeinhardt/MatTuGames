function gb=p_game_basis(n)
% P_GAMES_BASIS computes the game basis for n-players.
% Using Matlab's PCT.
%
% Usage: gb=p_game_basis(n)
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
%   05/21/2011        0.1 alpha        hme
%   06/20/2012        0.1 beta         hme
%   05/15/2014        0.5              hme
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


parfor k=1:N
   sSk=Subsets(k,n);
   gb(k,:)= ismembc(S,sSk);
end


%--------------------------------------
function sS=Subsets(S,n);

it=0:-1:1-n;
slcP=rem(floor(S(:)*pow2(it)),2)==0;

J=1:n;
sP=J(slcP);

S1=1:S; 

if (2^n-1)==S
  sS=S1;
else 
 for k=1:length(sP)
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
